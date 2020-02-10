"""
Helpers for handling data and metadata.
"""

import csv
import re
import logging
import hashlib
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
    """Get the raw data files for a single run, either from local disk or from URLs.

    runid: Illumina Run ID
    outdir: run directory to save fastq.gz to
    runs: dictionary of per-run metadata
    """
    LOGGER.info("get_data: runid %s", runid)
    LOGGER.info("get_data: outdir %s", outdir)
    rundir = Path("/seq/runs") / runid
    outfiles = {
        "R1": Path(outdir) / "Undetermined_S0_L001_R1_001.fastq.gz",
        "R2": Path(outdir) / "Undetermined_S0_L001_R2_001.fastq.gz",
        "I1": Path(outdir) / "Undetermined_S0_L001_I1_001.fastq.gz"}
    # Get data from raw run data
    if rundir.is_dir() and False:
        LOGGER.info("get_data: running bcl2fastq")
        shell(
            """
                bcl2fastq --create-fastq-for-index-reads \
                    -R /seq/runs/{runid} \
                    -o {outdir}
            """)
    # Get data from URLs
    else:
        url = runs[runid].get("URL")
        # If there's a unified URL, it's a zip file with R2/R2/I1 inside.
        if url:
            LOGGER.info("get_data: downloading single zip")
            shell("cd $(dirname {outdir}) && wget -q '{url}' && unzip \"$(basename '{url}')\"")
        # If there are separate URLs, they're the R1/R2/I1 fastq.gz themselves.
        else:
            urls = {key: runs[runid].get(key) for key in ("URLR1", "URLR2", "URLI1")}
            if urls["URLR1"] and urls["URLR2"] and urls["URLI1"]:
                LOGGER.info("get_data: downloading separate zips")
                urlr1 = urls["URLR1"]
                urlr2 = urls["URLR2"]
                urli1 = urls["URLI1"]
                outr1 = outfiles["R1"].name
                outr2 = outfiles["R2"].name
                outi1 = outfiles["I1"].name
                shell(
                    """
                        cd {outdir}
                        wget -q '{urlr1}' && mv $(basename '{urlr1}') {outr1}
                        wget -q '{urlr2}' && mv $(basename '{urlr2}') {outr2}
                        wget -q '{urli1}' && mv $(basename '{urli1}') {outi1}
                    """)
            else:
                raise ValueError("Need raw data or URLs for run %s" % runid)
    # Check files against expected checksums, if present.
    md5s = {key: runs[runid].get(key) for key in ("MD5R1", "MD5R2", "MD5I1")}
    if md5s["MD5R1"] and not md5s["MD5R1"] == md5(outfiles["R1"]):
        outfiles["R1"].rename(str(outfiles["R1"]) + ".failed")
        raise ValueError("MD5 mismatch on %s" % outfiles["R1"])
    if md5s["MD5R2"] and not md5s["MD5R2"] == md5(outfiles["R2"]):
        outfiles["R2"].rename(str(outfiles["R2"]) + ".failed")
        raise ValueError("MD5 mismatch on %s" % outfiles["R2"])
    if md5s["MD5I1"] and not md5s["MD5I1"] == md5(outfiles["I1"]):
        outfiles["I1"].rename(str(outfiles["I1"]) + ".failed")
        raise ValueError("MD5 mismatch on %s" % outfiles["I1"])

def md5(fp_in):
    """Get MD5 checksum of a file as a hex string."""
    # https://stackoverflow.com/a/3431838/4499968
    hash_md5 = hashlib.md5()
    with open(fp_in, "rb") as f_in:
        for chunk in iter(lambda: f_in.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
