"""
Helpers for handling data and metadata.
"""

import re
import csv
import logging
import hashlib
import gzip
import yaml
import urllib.request
from io import StringIO
from pathlib import Path
from Bio import SeqIO
from snakemake.shell import shell

LOGGER = logging.getLogger(__name__)

def download_metadata(path_yaml, path_out=None):
    """Download CSV from google sheets.

    YAML like this:
        metadata:
          slug: 2PACX-1vRsc...
          sheets:
            samples: 1692570416
            specimens: 1355081593
            ...
    Becomes:
         metadata/samples.csv
         metadata/specimens.csv
         metadata/...
    """
    with open(path_yaml) as f_in:
        info = yaml.safe_load(f_in)
    for key in info:
        slug = info[key]["slug"]
        sheets = info[key]["sheets"]
        for sheet_name, gid in sheets.items():
            url = f"https://docs.google.com/spreadsheets/d/e/{slug}/pub?gid={gid}&single=true&output=csv"
            path = Path(key)/(sheet_name + ".csv")
            if path_out is not None and str(path) == path_out:
                with urllib.request.urlopen(url) as f_in, open(path, "wt") as f_out:
                    reader = csv.reader(StringIO(f_in.read().decode("UTF8")))
                    writer = csv.writer(f_out, lineterminator="\n")
                    for row in reader:
                        writer.writerow(row)

class MetadataError(Exception):
    """Errors related to missing or invalid metadata."""

def load_runs(fp_in):
    """Load run metadata CSV.

    Output is a dictionary of run IDs to run attributes.
    """
    LOGGER.info("load_runs: fp_in %s", fp_in)
    runs = load_csv(fp_in, "Run")
    runs = {k: v for k, v in runs.items() if v["Skip"] != "TRUE"}
    return runs

def load_specimens(fp_in):
    """Load specimen metadata CSV.

    Output is a dictionary of specimen names to their attributes.
    """
    LOGGER.info("load_specimens: fp_in %s", fp_in)
    return load_csv(fp_in, "Specimen")

def load_samples(fp_in, specimens=None, runs=None):
    """Load sample metadata CSV and optionally link to barcode data.

    Output is a dictionary of sample names to sample attributes.  If
    dictionaries for specimens, runs, or sequences are given, those are joined
    to the sample attributes as nested dictionaries for each sample.
    """
    LOGGER.info("load_samples: fp_in %s", fp_in)
    LOGGER.info("load_samples: specimens %s...", str(specimens)[0:60])
    LOGGER.info("load_samples: runs %s...", str(runs)[0:60])
    samples = load_csv(fp_in, "Sample")
    samples = {k: v for k, v in samples.items() if v["Skip"] != "TRUE"}
    noblanks = {k: v for k, v in samples.items() if \
            v["BarcodeFwd"] and v["BarcodeRev"] and v["Run"]}
    if len(noblanks) < len(samples):
        LOGGER.warning("Blanks in BarcodeFwd/BarcodeRev/Run columns")
    unique_keys = {(row["BarcodeFwd"], row["BarcodeRev"], row["Run"]) for row in noblanks.values()}
    if len(unique_keys) != len(noblanks.values()):
        msg = "Duplicate entries in CSV %s for BarcodeFwd/BarcodeRev/Run combinations" % fp_in
        LOGGER.critical(msg)
        raise MetadataError(msg)
    for sample in samples.values():
        # Missing all three of these is a special case where we have the
        # preprocessed files but not the raw material
        if sample["BarcodeFwd"] or sample["BarcodeRev"] or sample["Run"]:
            if runs:
                _load_nested_items(sample, "Sample", runs, "Run")
        if specimens:
            _load_nested_items(sample, "Sample", specimens, "Specimen")
    return samples

def load_antibody_lineages(fp_in):
    """Load antibody lineage CSV."""
    LOGGER.info("load_antibody_lineages: fp_in %s", fp_in)
    lineages = load_csv(fp_in, "AntibodyLineage")
    return lineages

def load_antibody_isolates(fp_in, antibody_lineages=None):
    """Load antibody isolate CSV and optionally link to antibody lineage data.

    Output is a dictionary of isolate names to isolate attributes.  If a
    dictionary for lineages is given, entries will be joined as nested
    dictionaries.
    """
    LOGGER.info("load_antibody_isolates: fp_in %s", fp_in)
    LOGGER.info("load_antibody_isolates: antibody_lineages %s...", str(antibody_lineages)[0:60])
    isolates = load_csv(fp_in, "AntibodyIsolate")
    if antibody_lineages:
        for isolate in isolates.values():
            _load_nested_items(
                isolate, "AntibodyIsolate", antibody_lineages, "AntibodyLineage")
    return isolates

def _load_nested_items(entry, entrykey, others, key):
    """Add additional nested metdata to an entity from another dictionary."""
    other = others.get(entry[key])
    if other:
        entry[key + "Attrs"] = other.copy()
    else:
        LOGGER.error("Missing %s for %s %s", key, entrykey, entry[entrykey])

def load_csv(fp_in, key=None):
    """Generic CSV to dictionary.

    Assumes first row is header names.  If key is given, that column is used as
    the key for the top-level dictionary output.  If not, the first column is
    used.
    """
    LOGGER.debug("load_csv: fp_in %s", fp_in)
    warn = False
    with open(fp_in) as f_in:
        entries = {}
        reader = csv.DictReader(f_in)
        if not key:
            key = reader.fieldnames[0]
        for row in reader:
            row_key = row[key]
            if not row_key:
                warn = True
                continue
            if row_key in entries:
                msg = "Duplicate keys in CSV %s under column %s (value %s)" % (
                    fp_in, key, row_key)
                LOGGER.critical(msg)
                raise MetadataError(msg)
            entries[row_key] = row
    if warn:
        LOGGER.error(f"skipped rows with missing key value from {fp_in}")
    return entries

def get_data(runid, outdir, runs, runpath=None, threads=None):
    """Get the raw data files for a single run, either from local disk or from URLs.

    runid: Illumina Run ID
    outdir: run directory to save fastq.gz to
    runs: dictionary of per-run metadata
    runpath: path on local disk where run directories are found
    """
    LOGGER.info("get_data: runid %s", runid)
    LOGGER.info("get_data: outdir %s", outdir)
    if runpath:
        rundir = Path(runpath) / runid
    else:
        rundir = None
    outfiles = {
        "R1": Path(outdir) / "Undetermined_S0_L001_R1_001.fastq.gz",
        "R2": Path(outdir) / "Undetermined_S0_L001_R2_001.fastq.gz",
        "I1": Path(outdir) / "Undetermined_S0_L001_I1_001.fastq.gz"}
    # Get data from raw run data
    if rundir and rundir.is_dir():
        _get_data_bcl2fastq(runid, outdir, threads)
    # Get data from URLs
    else:
        _get_data_download(runid, outdir, runs, outfiles)
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

def _get_data_bcl2fastq(runid, outdir, threads):
    """Get data locally by running bcl2fastq on raw run directory."""
    LOGGER.info("get_data: running bcl2fastq")
    if threads:
        shell(
            """
                bcl2fastq --create-fastq-for-index-reads \
                    -r {threads} -d {threads} -p {threads} -w {threads} \
                    -R /seq/runs/{runid} \
                    -o {outdir}
            """.format(runid=runid, outdir=outdir, threads=threads))
    else:
        shell(
            """
                bcl2fastq --create-fastq-for-index-reads \
                    -R /seq/runs/{runid} \
                    -o {outdir}
            """.format(runid=runid, outdir=outdir))

def _get_data_download(runid, outdir, runs, outfiles):
    """Get data remotely by downloading zip files."""
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

def md5(fp_in):
    """Get MD5 checksum of a file as a hex string."""
    # https://stackoverflow.com/a/3431838/4499968
    hash_md5 = hashlib.md5()
    with open(fp_in, "rb") as f_in:
        for chunk in iter(lambda: f_in.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def get_samples_per_run(samples):
    """Make a dictionary of run names -> lists of sample names.

    Samples with no run specified are skipped.
    """
    samples_per_run = {}
    for sample_name, sample_attrs in samples.items():
        runid = sample_attrs["Run"]
        if not runid:
            continue
        if not runid in samples_per_run:
            samples_per_run[runid] = []
        samples_per_run[runid].append(sample_name)
    return samples_per_run

def amplicon_files(pattern, samples, cell_type_keep=None):
    """Make list of filenames for all specimen/chain/type combos of a given cell type.

    Give a pattern with templated fields for specimen, chain, and chain_type.
    Escape any others.  We're using the sample metadata for this, but it's not
    really sample-specific as we may have more than one sample for the same
    specimen/chain/chain type.  We'll deduplicate at the end in case.
    """
    LOGGER.info("amplicon_files: pattern: %s", pattern)
    LOGGER.info("amplicon_files: samples: %d entries", len(samples.keys()))
    LOGGER.info("amplicon_files: cell_type_keep: %s", cell_type_keep)
    target = []
    for samp_name, samp_items in samples.items():
        spec_name = samp_items["Specimen"]
        cell_type = samp_items["SpecimenAttrs"]["CellType"]
        chain = samp_items["Chain"]
        chain_type = samp_items["Type"]
        LOGGER.debug(
            "amplicon_files: process sample: %s " \
            "(specimen %s, cell type %s, chain %s, chain type %s)",
            samp_name, spec_name, cell_type, chain, chain_type)
        if not cell_type_keep or cell_type_keep in cell_type:
            text = pattern.format(chain=chain, chain_type=chain_type, specimen=spec_name)
            if text not in target:
                LOGGER.debug("amplicon_files: add target: %s", text)
                target.append(text)
    return target

def transpose_sample_md(samples, celltype):
    """Set up dictionaries of lists of sample metadata for a cell type.

    NOTE: I think this overlaps a lot with amplicon_files.  refactor?
    """
    sonar_md = {
        "chains": [sample["Chain"] for sample in samples.values()],
        "chaintypes": [sample["Type"] for sample in samples.values()],
        "celltypes": [sample["SpecimenAttrs"]["CellType"] for sample in samples.values()],
        "subjects": [sample["SpecimenAttrs"]["Subject"] for sample in samples.values()],
        "timepoints": [sample["SpecimenAttrs"]["Timepoint"] for sample in samples.values()],
        "specimens": [sample["Specimen"] for sample in samples.values()]
    }
    keeps = [celltype in entry for entry in sonar_md["celltypes"]]
    for key in sonar_md:
        sonar_md[key] = [val for i, val in enumerate(sonar_md[key]) if keeps[i]]
    return sonar_md
