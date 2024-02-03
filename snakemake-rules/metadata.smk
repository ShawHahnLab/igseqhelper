"""
Handlers for metadata.
"""

import csv
import logging
import yaml
import urllib.request
from io import StringIO
from pathlib import Path

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
        raise ValueError(msg)
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
                raise ValueError(msg)
            entries[row_key] = row
    if warn:
        LOGGER.error(f"skipped rows with missing key value from {fp_in}")
    return entries

rule all_get_metadata:
    input:
       specimens="metadata/specimens.csv",
       runs="metadata/runs.csv",
       samples="metadata/samples.csv",
       antibody_lineages="metadata/antibody_lineages.csv",
       antibody_isolates="metadata/antibody_isolates.csv"

def _setup_metadata(fp_specimens, fp_runs, fp_samples, fp_antibody_lineages, fp_antibody_isolates):
    global SPECIMENS, RUNS, SAMPLES, ANTIBODY_LINEAGES, ANTIBODY_ISOLATES
    SPECIMENS = load_specimens(fp_specimens)
    RUNS = load_runs(fp_runs)
    SAMPLES = load_samples(
        fp_samples,
        specimens=SPECIMENS,
        runs=RUNS)
    ANTIBODY_LINEAGES = load_antibody_lineages(fp_antibody_lineages)
    ANTIBODY_ISOLATES = load_antibody_isolates(fp_antibody_isolates, ANTIBODY_LINEAGES)

try:
    _setup_metadata(
        "metadata/specimens.csv",
        "metadata/runs.csv",
        "metadata/samples.csv",
        "metadata/antibody_lineages.csv",
        "metadata/antibody_isolates.csv")
except FileNotFoundError:
    print("Skipping metadata loading; be sure to run get_metadata rule.")
    SPECIMENS = {}
    RUNS = {}
    SAMPLES = {}
    ANTIBODY_LINEAGES = {}
    ANTIBODY_ISOLATES = {}

rule get_metadata:
    """Create CSV files from Google sheets based on metadata YAML."""
    output: "metadata/{sheet}.csv"
    input: "metadata.yml"
    run:
        igseqhelper.data.download_metadata(input[0], output[0])
