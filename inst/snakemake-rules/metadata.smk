"""
Handlers for metadata.
"""

import igseqhelper.data

rule all_get_metadata:
    input:
       specimens="metadata/specimens.csv",
       runs="metadata/runs.csv",
       samples="metadata/samples.csv",
       antibody_lineages="metadata/antibody_lineages.csv",
       antibody_isolates="metadata/antibody_isolates.csv"

def _setup_metadata(fp_specimens, fp_runs, fp_samples, fp_antibody_lineages, fp_antibody_isolates):
    global SPECIMENS, RUNS, SAMPLES, ANTIBODY_LINEAGES, ANTIBODY_ISOLATES
    SPECIMENS = igseqhelper.data.load_specimens(fp_specimens)
    RUNS = igseqhelper.data.load_runs(fp_runs)
    SAMPLES = igseqhelper.data.load_samples(
        fp_samples,
        specimens=SPECIMENS,
        runs=RUNS)
    ANTIBODY_LINEAGES = igseqhelper.data.load_antibody_lineages(fp_antibody_lineages)
    ANTIBODY_ISOLATES = igseqhelper.data.load_antibody_isolates(fp_antibody_isolates, ANTIBODY_LINEAGES)

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
