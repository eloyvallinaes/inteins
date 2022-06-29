"""
Scan sequences for prosite motifs
"""

import subprocess

SIGNATURES = ["TIGR01443.1", "TIGR01445.1"] #cterm, nterm


def scan(cogcode):
    """
    Scan sequences in fasta format with HMM of the intein signatures.
    Requieres hmmer installed or in PATH to use hmmsearch.
    """
    for signature in SIGNATURES:
        subprocess.run([
                "hmmsearch",
                # "-E",
                # "50",
                # "--domE",
                # "50",
                "--domtblout",
                f"scans/{signature}_{cogcode}.tab",
                f"hmm/{signature}.HMM",
                f"fasta/{cogcode}.fasta",
            ],
            stdout=subprocess.DEVNULL
        )


if __name__ == "__main__":
    for cogid in ["COG0209", "COG0305", "COG0417"]:
        scan(cogid)
