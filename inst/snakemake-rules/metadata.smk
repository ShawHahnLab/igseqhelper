"""
Handlers for metadata.
"""

import igseq.data

rule all_get_metadata:
    input:
       sequences="metadata/sequences.csv",
       specimens="metadata/specimens.csv",
       runs="metadata/runs.csv",
       samples="metadata/samples.csv",
       antibody_lineages="metadata/antibody_lineages.csv",
       antibody_isolates="metadata/antibody_isolates.csv"

def _setup_metadata(fp_sequences, fp_specimens, fp_runs, fp_samples, fp_antibody_lineages, fp_antibody_isolates):
    global SEQUENCES, SPECIMENS, RUNS, SAMPLES, ANTIBODY_LINEAGES, ANTIBODY_ISOLATES
    SEQUENCES = igseq.data.load_sequences(fp_sequences)
    SPECIMENS = igseq.data.load_specimens(fp_specimens)
    RUNS = igseq.data.load_runs(fp_runs)
    SAMPLES = igseq.data.load_samples(
        fp_samples,
        specimens=SPECIMENS,
        runs=RUNS,
        sequences=SEQUENCES)
    ANTIBODY_LINEAGES = igseq.data.load_antibody_lineages(fp_antibody_lineages)
    ANTIBODY_ISOLATES = igseq.data.load_antibody_isolates(fp_antibody_isolates, ANTIBODY_LINEAGES)

try:
    _setup_metadata(
        "metadata/sequences.csv",
        "metadata/specimens.csv",
        "metadata/runs.csv",
        "metadata/samples.csv",
        "metadata/antibody_lineages.csv",
        "metadata/antibody_isolates.csv")
except FileNotFoundError:
    print("Skipping metadata loading; be sure to run get_metadata rule.")
    SEQUENCES = {}
    SPECIMENS = {}
    RUNS = {}
    SAMPLES = {}
    ANTIBODY_LINEAGES = {}
    ANTIBODY_ISOLATES = {}

checkpoint get_metadata:
    """Create CSV files from Google sheets based on metadata YAML."""
    output:
       sequences="metadata/sequences.csv",
       specimens="metadata/specimens.csv",
       runs="metadata/runs.csv",
       samples="metadata/samples.csv",
       antibody_lineages="metadata/antibody_lineages.csv",
       antibody_isolates="metadata/antibody_isolates.csv"
    input: "metadata.yml"
    run:
        R("""
          devtools::load_all("igseq")
          update_metadata_via_yaml("{input}", ".")
          """)
        _setup_metadata(
            output.sequences,
            output.specimens,
            output.runs,
            output.samples,
            output.antibody_lineages,
            output.antibody_isolates)
