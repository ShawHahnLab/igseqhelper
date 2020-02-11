"""
Handlers for data/metadata
"""

from pathlib import Path
from igseq.data import (load_samples, load_specimens, load_runs, load_sequences, get_data)

def _setup_metadata(fp_samples, fp_specimens, fp_runs, fp_sequences):
    global SEQUENCES, SPECIMENS, SAMPLES, RUNS
    # key is sequence name, entries are dicts
    SEQUENCES = load_sequences(fp_sequences)
    # key is specimen name, entries are dicts
    SPECIMENS = load_specimens(fp_specimens)
    # key is run, entries are dicts
    RUNS = load_runs(fp_runs)
    # key is sample, entries are nested dictionaries of sample attributes and
    # other info above
    SAMPLES = load_samples(
        fp_samples,
        specimens=SPECIMENS,
        runs=RUNS,
        sequences=SEQUENCES)

try:
    _setup_metadata(
        "metadata/samples.csv",
        "metadata/specimens.csv",
        "metadata/runs.csv",
        "metadata/sequences.csv")
except FileNotFoundError:
    print("Skipping metadata loading; be sure to run get_metadata rule.")
    SEQUENCES = {}
    SAMPLES = {}
    RUNS = {}

SAMPLES_H = []
SAMPLES_L = []
ALLELES = {"heavy": {"V": "", "D": "", "J": ""}, "light": {"V": "", "J": ""}}

RAW = "Undetermined_S0_L001_{rp}_001.fastq.gz"


rule all_get_data:
    """Get data files for all sequencing runs.

    This lists the metadata checkpoint output as input so that we know what the
    sequencing runs are before running the rule.
    """
    input:
        files=expand("data/{run}/" + RAW, run = SAMPLES.keys(), rp = ["R1", "R2", "I1"]),
        metadata=lambda w: checkpoints.get_metadata.get().output

rule all_get_metadata:
    input:
       samples="metadata/samples.csv",
       specimens="metadata/specimens.csv",
       sequences="metadata/sequences.csv",
       runs="metadata/runs.csv"

rule get_data:
    """Get data files for a run.

    If raw data is available locally, run Illumina's bcl2fastq to create R1/R2/I1
    files.  (Add location of bcl2fatq to PATH if needed.)  If not, download from
    URLs listed in metadata.
    """
    output: expand("data/{run}/" + RAW, run = "{run}", rp = ["R1", "R2", "I1"])
    input: lambda w: checkpoints.get_metadata.get().output
    run: get_data(wildcards.run, Path(output[0]).parent, RUNS)

checkpoint get_metadata:
    """Create CSV files from Google sheets based on metadata YAML."""
    output:
       samples="metadata/samples.csv",
       specimens="metadata/specimens.csv",
       runs="metadata/runs.csv",
       sequences="metadata/sequences.csv"
    input: "metadata.yml"
    run:
        R("""
          devtools::load_all("igseq")
          update_metadata_via_yaml("{input}", ".")
          """)
        _setup_metadata(output.samples, output.specimens, output.runs, output.sequences)
