"""
Fetch sequences datasets from different sources
"""

import csv
import math
import requests
import pandas as pd
from Bio import Entrez
Entrez.email = "eloyvallina33@gmail.com"

class RefSeqFetch:
    taxkeys = ["kingdom", "phylum", "subphylum", "class", "order", "family", "genus"]
    def set_ids(self, ids):
        self.ids = ids

    def getdata(self):
        handle = self.fetch()
        rows = []
        for entry in Entrez.read(handle):
            row = dict(
                            zip(
                                RefSeqFetch.taxkeys,
                                entry['GBSeq_taxonomy'].split("; ")
                            )
                        )
            row["refseq_accno"] = entry['GBSeq_accession-version']
            row["sequence"] = entry['GBSeq_sequence']
            rows.append(row)
        return rows

    def post(self):
        handle = Entrez.epost("protein", id=",".join(self.ids))
        return Entrez.read(handle).values()

    def fetch(self):
        query_key, WebEnv = self.post()
        handle = Entrez.efetch(
                        db="protein",
                        query_key=query_key,
                        WebEnv=WebEnv,
                        retmode="xml"
                        )
        return handle


class COGFetch:
    URL = "https://www.ncbi.nlm.nih.gov/research/cog/api/cog/?cog={}&format=json&page={}"
    def __init__(self, cogid: str, pagefrom=1, pageto=2, sequencer=RefSeqFetch()):
        self.cogid = cogid
        self.pagefrom = pagefrom
        self.pageto = pageto
        self.sequencer = sequencer
        self.records = self.getrecords()
        self.sequences = self.getseqdata()

    def readpage(self, page):
        """
        Send a GET request to the COG database URL and returns the content as JSON.
        """
        r = requests.get(COGFetch.URL.format(self.cogid, page))
        if r.ok:
            return r.json()

    def readrange(self):
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

    def getrecords(self, filename=""):
        records = []
        for record in self.readrange():
            records.append(
                dict(
                    evalue= record["evalue"],
                    taxid=  record["organism"]["taxid"],
                    assembly_id= record["organism"]["assembly_id"],
                    refseq_accno = record["protein"]["name"],
                    membership= record["membership"]["membership_class"],
                    funcats = "".join([
                                        cat["name"]
                                        for cat in record['cog']['funcats']
                                    ]),
                    cogid = record['cog']['cogid'],
                )
            )
        return records

    def getseqdata(self):
        ids = [r["refseq_accno"] for r in self.records]
        self.sequencer.set_ids(ids)
        return self.sequencer.getdata()


    def writeout(self, filename=""):
        seqdf = pd.DataFrame(self.sequences)
        cogdf = pd.DataFrame(self.records)
        out = seqdf.merge(cogdf, on = "refseq_accno")
        if filename:
            out.to_csv(filename, index=False)
            return

        return out





if __name__ == '__main__':
    fetcher = COGFetch("COG0417", pagefrom=1, pageto=2)
    df = fetcher.writeout()
