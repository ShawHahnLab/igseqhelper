### Reporting

from igseq.reporting import counts_sample_summary, counts_run_summary


TARGET_REPORT_COUNTS = []
for runid in RUNS.keys():
    TARGET_REPORT_COUNTS.extend(expand(
        "counts/demux/{run}/{sample}.{rp}.fastq.gz.counts",
        sample=SAMPLES_PER_RUN[runid] + ["unassigned"],
        run=runid,
        rp=["R1", "R2", "I1"]))

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
    run: counts_sample_summary(input[0], output[0], SAMPLES, RUNS, SAMPLES_PER_RUN)

rule counts_run_summary:
    """A per-run summary of raw read counts."""
    output: "reporting/counts_by_run.csv"
    input: "reporting/counts_by_sample.csv"
    run: counts_run_summary(input[0], output[0])
