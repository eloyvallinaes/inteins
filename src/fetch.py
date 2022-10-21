"""
Classes for fetching sequence data from different sources.
"""

import os
import csv
import math
import requests
import pandas as pd
from Bio import Entrez
Entrez.email = "eloyvallina33@gmail.com"


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


class RefSeqFetch:
    """
    Given RefSeq protein ids, retrieve sequence and taxonomic information
    over the NCBI Entrez interface.
    """
    taxkeys = [
        "kingdom",
        "phylum",
        "subphylum",
        "class",
        "order",
        "family",
        "genus"
    ]

    def set_ids(self, ids):
        self.ids = ids

    def post(self):
        for chunk in chunks(self.ids, 200):
            handle = Entrez.epost("protein", id=",".join(chunk))
            yield Entrez.read(handle).values()

    def fetch(self):
        for query_key, WebEnv in self.post():
            handle = Entrez.efetch(
                            db="protein",
                            query_key=query_key,
                            WebEnv=WebEnv,
                            retmode="xml"
                            )
            yield handle

    def getdata(self):
        rows = []
        for handle in self.fetch():
            for entry in Entrez.read(handle):
                row = dict(
                            zip(
                                RefSeqFetch.taxkeys,
                                entry['GBSeq_taxonomy'].split("; ")
                            )
                        )
                row["refseq_accno"] = entry['GBSeq_accession-version']
                row["sequence"] = entry['GBSeq_sequence'].upper()
                rows.append(row)
        return rows


class COGFetch:
    URL = "https://www.ncbi.nlm.nih.gov/research/cog/api/cog/?cog={}&format=json&page={}"

    def __init__(self, cogid: str, pagefrom=1, pageto=2, sequencer=RefSeqFetch()):
        self.cogid = cogid
        self.pagefrom = pagefrom
        self.pageto = pageto
        self.sequencer = sequencer
        self.records = self.getrecords()
        self.sequences = self.getseqdata()

    def readpage(self, page) -> dict:
        """
        Send a GET request to the COG database URL and return the content as JSON.
        """
        r = requests.get(COGFetch.URL.format(self.cogid, page))
        if r.ok:
            return r.json()

    def readrange(self) -> list[dict]:
        """
        Visit the pages in range pagefrom->pageto and collect results in data attribute.
        """
        data = []
        i = self.pagefrom
        while True:
            content = self.readpage(i)
            data.extend(content['results'])
            i += 1
            if content['next'] == None or i > self.pageto:
                break
        return data

    def getrecords(self) -> list[dict]:
        """
        Slice the COG JSON data obtained from a range of pages and return
        a list of dictionaries with fields of interest.
        """

        records = []
        for record in self.readrange():
            records.append(
                dict(
                    evalue=record["evalue"],
                    taxid=record["organism"]["taxid"],
                    assembly_id=record["organism"]["assembly_id"],
                    refseq_accno=record["protein"]["name"],
                    membership=record["membership"]["membership_class"],
                    funcats="".join([
                                        cat["name"]
                                        for cat in record['cog']['funcats']
                                    ]),
                    cogid=record['cog']['cogid'],
                )
            )
        return records

    def getseqdata(self) -> list[dict]:
        """
        Infer sequence accesions from COG data and obtaine sequences and
        taxonomic categoricals
        """
        ids = [r["refseq_accno"] for r in self.records]
        self.sequencer.set_ids(ids)
        return self.sequencer.getdata()

    def writeout(self, dirname="cogs"):
        """
        Merge sequence and COG data and write dataframe to file.
        """
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        out = self.toDataFrame()
        out.to_csv(f"{dirname}/{self.cogid}.tsv", index=False, sep="\t")

    def writefasta(self, dirname="fasta"):
        """
        Write sequences in fasta format.
        """
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        with open(f"{dirname}/{self.cogid}.fasta", "w") as fastafile:
            for row in self.sequences:
                acc = row['refseq_accno']
                seq = row['sequence']

                fastafile.write(f">{acc}\n{seq}\n")

    def toDataFrame(self) -> pd.DataFrame:
        """
        Merge sequence and COG data and return a dataframe
        """
        seqdf = pd.DataFrame(self.sequences)
        cogdf = pd.DataFrame(self.records)
        return cogdf.merge(seqdf, on = "refseq_accno")


if __name__ == '__main__':
    for cogid in ["COG0305", "COG0209", "COG0417"]:
        fetcher = COGFetch(cogid, pagefrom=1, pageto=math.inf)
        fetcher.writeout()
        fetcher.writefasta()
