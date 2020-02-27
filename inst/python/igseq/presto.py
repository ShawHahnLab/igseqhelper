"""
Helpers for running pRESTO.

This stores the configuration settings for the different pRESTO scripts and
provides helper functions for preparing input files.
"""

import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from igseq.data import load_sequences

LOGGER = logging.getLogger(__name__)

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
        "mean_qual": 30,
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
    LOGGER.info("prep_primers_fwd: fp_csv_in: %s", fp_csv_in)
    LOGGER.info("prep_primers_fwd: fp_fwd_out: %s", fp_fwd_out)
    sequences = load_sequences(fp_csv_in)
    fwd = sequences["5PIIA"]["Seq"]
    LOGGER.debug("prep_primers_fwd: fwd primer seq: %s", fwd)
    fwd = SeqRecord(Seq(fwd), id="5PIIA", description="")
    SeqIO.write(fwd, fp_fwd_out, "fasta")
