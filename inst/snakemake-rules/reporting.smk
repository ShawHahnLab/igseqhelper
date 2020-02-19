### Reporting

import igseq.reporting
from igseq.data import amplicon_files

TARGET_REPORT_ALL = expand(
    "reporting/{thing}.csv",
    thing=["counts_by_sample", "counts_by_run", "counts_amplicon_summary"])

TARGET_REPORT_COUNTS = expand(
    outputs_per_run("counts/demux/{run}/{sample}.{{rp}}.fastq.gz.counts", SAMPLES),
    rp=["R1", "R2", "I1"]) + expand(
    "counts/demux/{run}/unassigned.{rp}.fastq.gz.counts", run=RUNS.keys(), rp=["R1", "R2", "I1"])

TARGET_AMPLICON_COUNTS = amplicon_files(
    "counts/presto/data/{chain}.{chain_type}/{specimen}.R1.fastq.counts", SAMPLES)

TARGET_ASSEMBLY_COUNTS = amplicon_files(
    "counts/presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq.counts", SAMPLES)

TARGET_PRESTO_QUAL_COUNTS = amplicon_files(
    "counts/presto/qual/{chain}.{chain_type}/{specimen}_quality-pass.fastq.counts", SAMPLES)

rule report_all:
    input: TARGET_REPORT_ALL

# Sample-based

rule counts_table:
    """Just a big ol' list of file paths and read counts."""
    output: "reporting/counts.csv"
    input: TARGET_REPORT_COUNTS
    shell:
        """
            echo "Filename,NumSequences" > {output}
            paste -d , <(echo "{input}" | tr ' ' '\\n') <(cat {input}) >> {output}
        """

rule counts_sample_summary:
    """A per-sample summary of raw read counts."""
    output: "reporting/counts_by_sample.csv"
    input: "reporting/counts.csv"
    run: igseq.reporting.counts_sample_summary(input[0], output[0], SAMPLES)

rule counts_run_summary:
    """A per-run summary of raw read counts."""
    output: "reporting/counts_by_run.csv"
    input: "reporting/counts_by_sample.csv"
    run: igseq.reporting.counts_run_summary(input[0], output[0])

# Specimen+chain -based

rule counts_presto_amplicon_summary:
    """A per-amplicon summary of read counts."""
    output: "reporting/counts_amplicon_summary.csv"
    input: TARGET_AMPLICON_COUNTS
    run: igseq.reporting.counts_specimen_summary(input, output[0], SPECIMENS)

rule counts_presto_assembly_summary:
    """A per-specimen summary of paired read counts."""
    output: "reporting/counts_assembly_summary.csv"
    input: TARGET_ASSEMBLY_COUNTS
    run: igseq.reporting.counts_assembly_summary(input, output[0], SPECIMENS)

rule counts_presto_qual_summary:
    """A per-specimen summary of paired read counts."""
    output: "reporting/counts_presto_qual_summary.csv"
    input: TARGET_PRESTO_QUAL_COUNTS
    run: igseq.reporting.counts_presto_qual_summary(input, output[0], SPECIMENS)

# TODO next: also tally after pRESTO's QC and primer checking.  Maybe do a
# summary table in the report of counts following assembly, qc, and primer
# checking showing attrition across steps.
