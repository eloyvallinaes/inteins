"""
Parsers for different file formats used across this project.
"""

import re
import csv
import json
import copy
import pandas as pd
from pathlib import Path
from Bio import SeqIO, SearchIO
from typing import Iterator


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
    def __init__(self, filename):
        self.records = self.parse(filename)

    def parse(self, filename):
        with open(filename, "r") as infile:
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
    def __init__(self, filename, subset={}):
        self.subset = subset
        self.sequences, self.taxids, self.limits = self.parse(filename)
        self.annots = self.extract_annotations()
        self.hosts = self.extract_hosts()

    def parse(self, filename):
        """
        Load FASTA file to a dict as id:seq pairs. Optionally, limit parsing to
        IDs in subset.

        :param filename: output file path to write to.
        :type filename: str
        :param subset: a set of IDs to retrieve from FASTA file
        :type subset: set
        """
        sequences = dict()
        taxids = dict()
        limits = dict()
        with open(filename) as fastafile:
            for row in SeqIO.parse(fastafile, "fasta"):
                acc = row.id.split("|")[0].strip(">")
                taxid = row.id.split("|")[1]
                limstr = row.id.split("|")[2]
                if len(self.subset) == 0 or acc in self.subset:
                    sequences[acc] = row.seq.upper()
                    taxids[acc] = taxid
                    limits[acc] = Fasta2Dict.parse_limits(limstr)

        return sequences, taxids, limits

    def write_sequences(self, filename: str, use_taxids=False, subset={}):
        """
        Write object's sequences back to FASTA format. Optionally, limit
        output to a selection of IDs and/or write taxids as sequence headers.

        :param filename: output file path to write to.
        :type filename: str
        :param subset: a set of IDs to write to FASTA
        :type subset: set
        :param use_taxids: use taxids as sequence headers in FASTA file.
        :type use_taxids: bool

        """
        with open(filename, "w") as fastafile:
            for acc, seq in self.sequences.items():
                if acc not in subset and len(subset) > 0:
                    continue
                if use_taxids:
                    fastafile.write(
                            f">{self.taxids[acc]}\n{seq}\n"
                    )
                else:
                    fastafile.write(
                            f">{acc}\n{seq}\n"
                    )

    def write_segments(
                    self, filename: str, entity: str, use_taxids=False,
                    subset={}
                ):
        """
        Write FASTA file of annotations or hosts. Optionally replace accessions
        with taxids as sequence identifiers.

        :param filename: output destination file
        :type filename: str or pathlike
        :param entity: either 'annotations' or 'hosts'
        :type entity: str
        :param use_taxids: whether to replace accessions with taxids. Two
        sequences migh be given the same taxid, which can cause errors
        downstream.
        :type use_taxids: bool

        """
        if entity not in ["annots", "hosts"]:
            raise ValueError(
                f"{entity} not understood; must be 'inteins' or 'hosts'"
            )

        segments = self.annots if entity == "annots" else self.hosts
        with open(filename, "w") as fastafile:
            for acc, records in segments.items():
                if acc not in subset and len(subset) > 0:
                    continue
                for i, part in enumerate(records):
                    seq = part['seq']
                    if use_taxids:
                        taxid = self.taxids[acc]
                        fastafile.write(f">{taxid}|{entity}_{i}\n{seq}\n")

                    else:
                        fastafile.write(f">{acc}|{entity}_{i}\n{seq}\n")

    @staticmethod
    def parse_limits(limits: str):
        """
        Extract start and end positions from sequence headers regular language.

        :param limits: the string contaning the limits specified as integeres
        separated by '...', eg. IPR036844(3...109;139...177).
        :type limits: str
        :return: parsed limits, eg.
        [{'start': 3, 'end': 109}, {'start': 139, 'end': 177}]
        :rtype: list[dict], eg.


        """
        pattern = re.compile(r"([0-9]+)\.\.\.([0-9]+)")
        return [
            {'start': int(start), 'end': int(end)}
            for start, end in pattern.findall(limits)
        ]

    def extract_annotations(self):
        """
        Work out annotation segments from limits.
        """
        records = {}
        for acc, limEntries in self.limits.items():
            inteins = []
            for lim in limEntries:
                start = lim['start']
                end = lim['end']
                inteins.append({
                        'start': start,
                        'end': end,
                        'span': end-start,
                        'seq': self.sequences[acc][start:end]
                })

            records[acc] = inteins

        return records

    def extract_hosts(self):
        """
        Work out host protein spliced product from annotation limits.
        """
        hosts = dict()
        for acc, seq in self.sequences.items():
            indeces = [0]
            for record in self.limits[acc]:
                indeces += [record['start']] + [record['end']]
            indeces += [-1]

            host_seq = ''
            for e in range(0, len(indeces), 2):
                host_seq += self.sequences[acc][indeces[e]:indeces[e + 1]]

            hosts[acc] = [{
                "seq": host_seq, 'span': len(host_seq)
            }]

        return hosts

    def __sub__(self, other):
        """
        Implement Fasta2Dict subtraction as obtaing the difference between
        a long and short annotation, eg., to generate mininteins from an
        object of intein annotations minus an object of endonuclease
        annotations.

        """
        new = copy.copy(self)  # shallow copy is enough?
        selfKeys = self.annots.keys()
        otherKeys = other.annots.keys()
        interKeys = set(selfKeys).intersection(otherKeys)
        differences = dict()
        for acc in interKeys:
            differences[acc] = []
            for diff in Fasta2Dict.sequence_difference(self, other, acc):
                differences[acc].append({"seq": diff, "span": len(diff)})

        new.annots = differences
        return new

    def __add__(self, other):
        """
        Implement subtraction as obtaing the combination of two annotations
        such that two discontinous segments may be joined. Since we are adding
        sequences, this operation is generally non-commutative: a + b =/= b + a.

        """
        new = copy.copy(self)
        selfKeys = self.annots.keys()
        otherKeys = other.annots.keys()
        interKeys = set(selfKeys).intersection(otherKeys)
        joined = dict()
        for acc in interKeys:
            joined[acc] = []
            for cont in Fasta2Dict.sequence_join(self, other, acc):
                joined[acc].append({"seq": cont, "span": len(cont)})

        new.annots = joined
        return new

    @staticmethod
    def sequence_difference(outer, inner, acc) -> Iterator[str]:
        seq = outer.sequences[acc]
        for outrec, inrec in zip(outer.annots[acc], inner.annots[acc]):
            # enforce outer surrounds inner annotation
            if (
                outrec['start'] <= inrec['start'] and
                outrec['end'] >= inrec['end']
            ):
                diff = (
                    seq[outrec['start']:inrec['start']] +
                    seq[inrec['end']:outrec['end']]
                )
                yield diff
            else:
                yield ''

    @staticmethod
    def sequence_join(first, second, acc) -> Iterator[str]:
        for rec1, rec2 in zip(first.annots[acc], second.annots[acc]):
            yield rec1['seq'] + rec2['seq']


class Segments:
    @staticmethod
    def extractHost(sequence, inteins):
        indeces = [0]
        for i in inteins:
            indeces += [i['start']] + [i['end']]
        indeces += [-1]

        host = ''
        for e in range(0, len(indeces), 2):
            host += sequence[indeces[e]:indeces[e + 1]]

        return [{"seq": host}]


class InterproSegments(Segments):
    def __init__(self, filename):
        self.name = Path(filename).stem
        self.segments = self.parse(filename)

    def to_fasta(self, entity, outfilename):
        if entity not in ["intein", "host"]:
            raise ValueError(
                "{entity} not understood; must be 'inteins' or 'host'"
            )
        with open(outfilename, "w") as fastafile:
            for acc, values in self.segments.items():
                for i, part in enumerate(values[entity]):
                    sequence = part['seq']
                    fastafile.write(f">{acc}|{entity}_{i}\n{sequence}\n")

    @staticmethod
    def parse(filename):
        records = {}
        with open(filename) as fastafile:
            for row in SeqIO.parse(fastafile, "fasta"):
                acc, limits, name = row.id.split("|")
                seq = str(row.seq.upper())
                inteins = InterproSegments.extractInteins(seq, limits)
                host = Segments.extractHost(seq, inteins)
                records[acc] = {
                    'inteins': inteins,
                    'host': host,
                    }
        return records

    @staticmethod
    def parselimits(limits: str):
        """
        Extract limits information from from string regular language.
        """
        pattern = re.compile(r"([0-9]+)\.\.\.([0-9]+)")
        for start, end in pattern.findall(limits):
            yield int(start), int(end)

    @staticmethod
    def extractInteins(sequence, limits: str):
        """
        Create list of intein segments based on limits.
        """
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


class COGSegments(Segments):
    @staticmethod
    def parse(seqfile, ntermfile, ctermfile):
        seqmap = COGSegments.get_seqmap(seqfile)
        segments = COGSegments.extract_segments(ntermfile, ctermfile)\
                              .to_dict(orient="records")
        data = {}
        for acc, sequence in seqmap.items():
            segments = [
                segment
                for segment in segments
                if segment['refseq_accno'] == acc
            ]
            inteins = COGSegments.extract_inteins(sequence, segments)
            hostseq = Segments.extract_host(sequence, inteins)
            data[acc] = {'inteins': inteins, 'host': hostseq}
        return data

    @staticmethod
    def extract_inteins(sequence, segments):
        inteins = []
        for segment in segments:
            start = segment['start']
            end = segment['end']
            inteins.append({
                'start': start,
                'end': end,
                'span': end - start,
                'seq': sequence[start:end]
            })
        return inteins

    @staticmethod
    def get_seqmap(seqfile):
        sequences = pd.read_csv(seqfile, sep="\t").to_dict(orient="records")
        return {row['refseq_accno']: row['sequence'] for row in sequences}

    @staticmethod
    def extract_segments(ntermfile, ctermfile):
        xcolumns = ["end", "motifname", "motifacc"]
        ycolumns = ["start", "motifname", "motifacc"]
        nterm = pd.read_csv(ntermfile).drop(xcolumns, axis="columns")
        cterm = pd.read_csv(ctermfile).drop(ycolumns, axis="columns")
        ss = nterm.merge(cterm, on="refseq_accno")
        ss['span'] = ss.end - ss.start
        # select index combinations giving positive, shortest span for each acc
        # and output sorted by starting positions -- proper sequence slicing
        return ss[ss.span > 0].sort_values("span")\
                              .drop_duplicates(
                                ["refseq_accno", "start"],
                                keep="first")\
                              .sort_values("start")


class CDHit:
    """
    Interpret CD-Hit clstr file and extract:
    1. Cluster compositions
    2. Cluster representatives

    Interface with Fasta2Dict to write FASTA subsets for clusters.
    """
    accPattern = r">([A-Z]{1}[A-z0-9]{4,})"
    repPattern = r"[0-9]{1,}\s*[0-9]{1,}aa, >([A-Z]{1}[A-z0-9]{4,})\.\.\. *"

    def __init__(self, filename):
        self.clusters = self.parse(filename)  # dict
        self.ncluster = len(self.clusters)
        self.reps = self.get_representatives(filename)

    @classmethod
    def parse(cls, filename):
        result = []
        with open(filename, "r") as clstrfile:
            while True:
                line = clstrfile.readline()
                if not line:
                    break
                elif line.startswith('#'):
                    continue
                elif line.startswith(">"):
                    members = []
                    while True:
                        record = clstrfile.readline()
                        if not record:
                            break
                        elif record.startswith(">"):
                            result.append(members)
                            members = []
                        else:
                            members.append(
                                re.search(cls.accPattern, record).group(1)
                            )

            result.append(members)
        return result

    @classmethod
    def get_representatives(cls, filename):
        reps = []
        with open(filename, "r") as clstrfile:
            while True:
                line = clstrfile.readline()
                if not line:
                    break
                if "*" in line:
                    reps.append(
                        re.match(cls.repPattern, line).group(1)
                    )
        return reps

    def rep2members(self):
        """
        A dictionary mapping cluster representatives to cluster members,
        including itself.
        """
        return {r: m for r, m in zip(self.reps, self.clusters)}

    def cluster2Fasta(self, fastafile, clusteri, outname, replacement=dict()):
        """
        Write FASTA file containing all sequences in a cluster.
        Optionally, replace accessions in headers with any other key, such as
        taxids, through the Fasta2Dict interface.

        :param fastafile: a file containing the sequences from which clusters
        where obtained.
        :type fastafile: str
        :param clusteri: the cluster index to write from self.clusters
        :type clusteri: int
        :param outname: output FASTA filename
        :type outname: str
        :param replacement: a dictionary mapping accessions to alternative keys.
        See Fasta2Dict.
        :type replacement:

        """
        seqs = Fasta2Dict(fastafile, self.clusters[clusteri])
        if replacement:
            seqs.replaceKeys(replacement)

        seqs.write(outname)

    def reps2Fasta(self, fastafile, outname, replacement=dict()):
        """
        Write FASTA file containing representative sequences.
        Optionally, replace accessions in headers with any other key, such as
        taxids, through the Fasta2Dict interface.

        :param fastafile: a file containing the sequences from which clusters
        where obtained.
        :type fastafile: str
        :param outname: output FASTA filename
        :type outname: str
        :param replacement: a dictionary mapping accessions to alternative keys.
        See Fasta2Dict.
        :type replacement:
        """
        seqs = Fasta2Dict(fastafile, self.reps)
        if replacement:
            seqs.replaceKeys(replacement)

        seqs.write(outname)

    def membership(self):
        return {
            m: i
            for i, members in enumerate(self.clusters)
            for m in members
        }


def split_fasta_file(fastafile):
    """
    Break sequence collecition in FASTA format into single sequence files.

    """
    for records in SeqIO.parse(fastafile, 'fasta'):
        with open(f"{records.id}.fasta", "w") as seqfile:
            seqfile.write(f">{records.id}\n{records.seq}")


if __name__ == '__main__':
    fasta = Path("fasta/interpro")
    hint = Fasta2Dict(fasta / 'IPR036844.fasta')
    endo = Fasta2Dict(fasta / 'IPR027434.fasta')
