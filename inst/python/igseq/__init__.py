#!/usr/bin/env python

"""
Immunoglobulin Sequence Data Analysis

Various Python helpers for our analysis here, particularly to use in Snakemake
rules.  See also the igseq R package.
"""

from pathlib import Path
import gzip
import re
import csv
import sys
import logging
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio import SeqIO

def __find_r_pkg():
    LOGGER.debug("finding R package")
    igseq_inst = Path(__file__).parent.parent.parent
    if igseq_inst.name == "inst":
        igseq = igseq_inst.parent
    else:
        igseq_inst = igseq_inst.parent
        igseq = igseq_inst
    igseq_inst = igseq_inst.resolve()
    igseq = igseq.resolve()
    LOGGER.debug("finding R package: inst: %s", str(igseq_inst))
    LOGGER.debug("finding R package: pkg: %s", str(igseq))
    return igseq_inst, igseq

LOGGER = logging.getLogger(__name__)
#logging.basicConfig(level=0)

R_PKG_INST, R_PKG_PATH = __find_r_pkg()

def load_sequences(fp_in):
    """Load sequence metadata CSV.

    Output is a dictionary of sequence names to sequences and their attributes.
    """
    with open(fp_in) as f_in:
        reader = csv.DictReader(f_in)
        seqs = {row["Name"]: row for row in reader}
    return seqs

def load_samples(fp_in, sequences=None):
    """Load sample metadata CSV and optionally link to barcode data.

    Output is a dict where key is run, entries are lists of sample dicts.
    """
    samples = {}
    with open(fp_in) as f_in:
        reader = csv.DictReader(f_in)
        for row in reader:
            if not row["Run"] in samples:
                samples[row["Run"]] = []
            samples[row["Run"]].append(row)
    if sequences:
        for run in samples:
            samps = samples[run]
            for samp in samps:
                bc_fwd = sequences.get(samp["BarcodeFwd"])
                bc_rev = sequences.get(samp["BarcodeRev"])
                samp["BarcodeFwdSeq"] = re.sub(" ", "", bc_fwd)
                samp["BarcodeRevSeq"] = re.sub(" ", "", bc_rev)
                samp["BarcodePair"] = "%s %s" % (samp["BarcodeFwdSeq"], samp["BarcodeRevSeq"])

    return samples

def load_runs(fp_in):
    """Load run metadata CSV."""
    with open(fp_in) as f_in:
        runs = {}
        reader = csv.DictReader(f_in)
        for row in reader:
            runs[row["Run"]] = row
    return runs

def trim_alignment_to(fasta_fp_in, fasta_fp_out, seq_id):
    """Trim all sequences in a FASTA to match a given sequence ID.

    This removes start and end regions of other sequences so they line up with
    one particular sequence in the alignment.
    """
    for record in SeqIO.parse(fasta_fp_in, "fasta"):
        if record.id == seq_id:
            match = re.match("^-*(.*[^-])-*$", str(record.seq))
            start = match.start(1)
            end = match.end(1)
            break
    else:
        raise Exception("Seq ID %s not found" % seq_id)
    with open(fasta_fp_out, "w") as f_out:
        for record in SeqIO.parse(fasta_fp_in, "fasta"):
            SeqIO.write(record[start:end], f_out, "fasta")
