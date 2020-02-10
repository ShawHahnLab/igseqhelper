"""
Helpers for runnin pRESTO.

This stores the configuration settings for the different pRESTO scripts and
provides helper functions for preparing input files.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from igseq import load_sequences

PRESTO_OPTS = {
    "assembly": {
        # The format of the sequence identifier which defines
        # shared coordinate information across paired ends.
        # (default: presto)
        "coord": "illumina",
        # Specify which read to reverse complement before
        # stitching. (default: tail)
        "rc": "tail",
        # Minimum sequence length to scan for overlap in de novo assembly.
        # (default: 8)
        # This sounds a bit low for our case, but when I manually check some of
        # these (in the 15 - 25 nt overlap range) they do look correct, so I'll
        # leave it alone for now.
        "minlen": 8,
        # Maximum sequence length to scan for overlap in de novo assembly.
        # (default: 1000)
        # The shortest expected sequence should dictate this.  With 2x309 we
        # could never even consider more than 309 overlap, and we don't expect
        # anything (heavy or light) shorter than around 300 nt anyway.
        # Also leaving this at the default for the moment.
        "maxlen": 1000
    },
    "qc": {
        "mean_qual": 20,
        "fwd_start": 0,
        "fwd_mode": "cut",
        "fwd_pf": "VPRIMER",
        "rev_start": 0,
        "rev_mode": "cut",
        "rev_pf": "CPRIMER"
    },
    "collapse": {
        "uf": "CPRIMER",
        "cf": "VPRIMER",
        "f": "DUPCOUNT",
        "num": 2
    }
}

def prep_primers_fwd(fp_csv_in, fp_fwd_out):
    """Take our primer CSV and create fwd FASTA for pRESTO."""
    sequences = load_sequences(fp_csv_in)
    fwd = SeqRecord(Seq(sequences["5PIIA"]), id="5PIIA", description="")
    SeqIO.write(fwd, fp_fwd_out, "fasta")
