"""
Handlers for data/metadata
"""

from pathlib import Path
import igseq.data

def _setup_metadata(fp_samples, fp_specimens, fp_runs, fp_sequences, fp_antibody_lineages, fp_antibody_isolates):
    global SEQUENCES, SPECIMENS, RUNS, SAMPLES, ANTIBODY_LINEAGES, ANTIBODY_ISOLATES
    # key is sequence name, entries are dicts
    SEQUENCES = igseq.data.load_sequences(fp_sequences)
    # key is specimen name, entries are dicts
    SPECIMENS = igseq.data.load_specimens(fp_specimens)
    # key is run, entries are dicts
    RUNS = igseq.data.load_runs(fp_runs)
    # key is sample, entries are nested dictionaries of sample attributes and
    # other info above
    SAMPLES = igseq.data.load_samples(
        fp_samples,
        specimens=SPECIMENS,
        runs=RUNS,
        sequences=SEQUENCES)
    ANTIBODY_LINEAGES = igseq.data.load_antibody_lineages(fp_antibody_lineages)
    ANTIBODY_ISOLATES = igseq.data.load_antibody_isolates(fp_antibody_isolates, ANTIBODY_LINEAGES)

try:
    _setup_metadata(
        "metadata/samples.csv",
        "metadata/specimens.csv",
        "metadata/runs.csv",
        "metadata/sequences.csv",
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

RAW = "Undetermined_S0_L001_{rp}_001.fastq.gz"

def outputs_per_run(pattern, samples):
    """Helper to produce output filenames specific to run/sample combos.

    Give a string pattern with {run} and {sample} and the dictionary of sample
    information and the right run/sample combinations will be filled in.  Be
    sure to escape any other template variables like {{rp}}.
    """
    samples_per_run = igseq.data.get_samples_per_run(samples)
    target = []
    for runid in samples_per_run.keys():
        target.extend(expand(
            pattern,
            sample=samples_per_run[runid],
            run=runid))
    return target

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
    run: igseq.data.get_data(wildcards.run, Path(output[0]).parent, RUNS)

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
