import subprocess

from Bio import Entrez
Entrez.email = "eloyvallina33@gmail.com"

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def directFetch(filename):
    post = subprocess.Popen(
                [
                    "epost",
                    "-db",
                    "protein",
                    "-input",
                    filename,
                ],
                stdout=subprocess.PIPE
            )
    fetch = subprocess.Popen(
                [
                    "efetch",
                    "-db",
                    "protein",
                    "-format",
                    "fasta",
                    "-mode",
                    "text",
                ],
                stdin=post.stdout,
                stdout=subprocess.PIPE
            )
    post.stdout.close()
    out = fetch.communicate()[0]
    return out

class RefSeqFetch:
    taxkeys = ["kingdom", "phylum", "subphylum", "class", "order", "family", "genus"]
    def __init__(self, filename):
        self.filename = filename
        self.ids = [acc.strip() for acc in open(filename, "r").readlines()]

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



if __name__ == "__main__":
    cogid = "COG0209"
    for attempt in range(1):
        idfile = f"idlists/{cogid}_ids.txt"
        oricounts = len(open(idfile, "r").readlines())
        # test bio.Entrez
        fetcher = RefSeqFetch(idfile)
        bioseqs = fetcher.getdata()
        # counts = len(bioseqs)

        # test eutils subprocess
        directseqs = directFetch(idfile)
        # directcounts = directseqs.decode().count(">")

        # print(f"Attempt {attempt}: {counts} sequences for {cogid} (n = {oricounts})")
