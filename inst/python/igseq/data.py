"""
Helpers for handling data and metadata.
"""

import csv
import re
import logging
from pathlib import Path
from snakemake.shell import shell

LOGGER = logging.getLogger(__name__)

def load_sequences(fp_in):
    """Load sequence metadata CSV.

    Output is a dictionary of sequence names to sequences and their attributes.
    """
    LOGGER.info("load_sequences: fp_in %s", fp_in)
    with open(fp_in) as f_in:
        reader = csv.DictReader(f_in)
        seqs = {row["Name"]: row for row in reader}
    return seqs

def load_samples(fp_in, sequences=None):
    """Load sample metadata CSV and optionally link to barcode data.

    Output is a dict where key is run, entries are lists of sample dicts.
    """
    LOGGER.info("load_samples: fp_in %s", fp_in)
    LOGGER.info("load_samples: sequences %s...", str(sequences)[0:60])
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
                samp["BarcodeFwdSeq"] = re.sub(" ", "", bc_fwd["Seq"])
                samp["BarcodeRevSeq"] = re.sub(" ", "", bc_rev["Seq"])
                samp["BarcodePair"] = "%s %s" % (samp["BarcodeFwdSeq"], samp["BarcodeRevSeq"])

    return samples

def load_runs(fp_in):
    """Load run metadata CSV."""
    LOGGER.info("load_runs: fp_in %s", fp_in)
    with open(fp_in) as f_in:
        runs = {}
        reader = csv.DictReader(f_in)
        for row in reader:
            runs[row["Run"]] = row
    return runs

def get_data(runid, outdir, runs):
    """Get the raw data files for a single run, either from local disk or from URLs."""
    LOGGER.info("get_data: runid %s", runid)
    LOGGER.info("get_data: outdir %s", outdir)
    rundir = Path("/seq/runs") / runid
    if rundir.is_dir():
        LOGGER.info("get_data: running bcl2fastq")
        shell(
            """
                bcl2fastq --create-fastq-for-index-reads \
                    -R /seq/runs/%s \
                    -o %s
            """ % (runid, outdir))
    else:
        url = runs[runid].get("URL")
        if url:
            LOGGER.info("get_data: downloading single zip")
            shell("cd data/{wildcards.run} && wget '%s' && unzip {wildcards.run}.zip" % url)
        else:
            url_r1 = runs[runid].get("URLR1")
            url_r2 = runs[runid].get("URLR2")
            url_i1 = runs[runid].get("URLI1")
            if url_r1 and url_r2 and url_i1:
                LOGGER.info("get_data: downloading separate zips")
                shell(
                    """
                        cd data/{wildcards.run}
                        wget '%s' && unzip %s_R1.zip
                        wget '%s' && unzip %s_R2.zip
                        wget '%s' && unzip %s_I1.zip
                    """ % (url_r1, runid, url_r2, runid, url_i1, runid))
            else:
                raise ValueError("Need raw data or URLs for run %s" % runid)
