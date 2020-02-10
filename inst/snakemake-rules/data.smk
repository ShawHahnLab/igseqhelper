"""
Handlers for data/metadata
"""

from pathlib import Path
from igseq.data import (load_sequences, load_samples, load_runs, get_data)

def _setup_metadata(fp_primers, fp_samples, fp_runs):
    global SEQUENCES, SAMPLES, RUNS, SAMPLES_ALL
    # key is sequence name, entries are dicts
    SEQUENCES = load_sequences("metadata/sequences.csv")
    # key is run, entries are lists of sample dicts
    SAMPLES = load_samples("metadata/samples.csv", SEQUENCES)
    # key is run, entries are dicts
    RUNS = load_runs("metadata/runs.csv")
    # A set of all unique samples
    SAMPLES_ALL = set()
    for samps in SAMPLES.values():
        SAMPLES_ALL = SAMPLES_ALL | set([s["Sample"] for s in samps])

try:
    _setup_metadata("metadata/sequences.csv", "metadata/samples.csv", "metadata/runs.csv")
except FileNotFoundError:
    print("Skipping metadata loading; be sure to run get_metadata rule.")
    SEQUENCES = {}
    SAMPLES = {}
    RUNS = {}
    SAMPLES_ALL = set()

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
       sequences="metadata/sequences.csv",
       runs="metadata/runs.csv"
    input: "metadata.yml"
    run:
        R("""
          devtools::load_all("igseq")
          update_metadata_via_yaml("{input}", ".")
          """)
        _setup_metadata(output.sequences, output.samples, output.runs)
