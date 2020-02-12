"""
Various QC checks.
"""

from igseq.qc import make_qualtrim_csv

TARGET_QUALTRIM_GRID = expand("qc/{run}/qualtrim.{sample}.{rp}.csv",
    sample = SAMPLES.keys(),
    run = RUNS.keys(),
    rp = ["R1", "R2"])

rule all_qualtrim_grid:
    input: TARGET_QUALTRIM_GRID

rule qualtrim_grid:
    """Make a CSV table summarizing cutadapt trim cutoffs vs output length."""
    output: "qc/{run}/qualtrim.{sample}.{rp}.csv"
    input: "demux/{run}/{sample}.{rp}.fastq.gz"
    run: make_qualtrim_csv(input[0], output[0])