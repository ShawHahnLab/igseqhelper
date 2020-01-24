"""
Helpers for runnin pRESTO.

This stores the configuration settings for the different pRESTO scripts and
provides helper functions for preparing input files.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from igseq import load_primers

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
    # Load in the primers, keeping only the forward cases, and trimming to the
    # portion after the barcode.
    fmtfwd = lambda s: s.lstrip("N")[8:]
    primers = load_primers(fp_csv_in)
    primers = {key: fmtfwd(primers[key]) for key in primers if "P5" in key}
    # If the're all the same, as I expect, just store one.
    if len(set(primers.values())) == 1:
        primers = {"Fwd": list(primers.values())[0]}
    mkseq = lambda key: SeqRecord(Seq(primers[key]), id=key, description="")
    fwd = [mkseq(key) for key in primers]
    SeqIO.write(fwd, fp_fwd_out, "fasta")

def prep_primers_full(fp_csv_in, fp_fwd_out, fp_rev_out):
    """Take our primer CSV and create fwd/rev FASTA for pRESTO.

    I don't think we actually have anything to use for the reverse case because
    of how our libraries were prepared.  See prep_primers_fwd instead.
    """
    primers = load_primers(fp_csv_in)
    fwd_keys = [key for key in primers if "P5_" in key]
    fmtfwd = lambda s: s.lstrip("N")[8:]
    fmtrev = lambda s: s
    mkseqfwd = lambda key: SeqRecord(Seq(fmtfwd(primers[key])), id=key, description="")
    mkseqrev = lambda key: SeqRecord(Seq(fmtrev(primers[key])), id=key, description="")
    fwd = [mkseqfwd(key) for key in fwd_keys]
    SeqIO.write(fwd, fp_fwd_out, "fasta")
    rev_keys = [key for key in primers if "P7_" in key]
    rev = [mkseqrev(key) for key in rev_keys]
    SeqIO.write(rev, fp_rev_out, "fasta")
