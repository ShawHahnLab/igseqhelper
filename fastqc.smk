rule fastqc_merge:
    input: expand("analysis/fastqc/merge/{run}/{sample}_fastqc.html", zip, run=[attrs["Run"] for attrs in SAMPLES.values() if attrs["Run"]], sample=[attrs["Sample"] for attrs in SAMPLES.values() if attrs["Run"]])

rule fastqc_reads:
    input: expand("analysis/fastqc/reads/{run}/Undetermined_S0_L001_{rp}_001_fastqc.html", run=RUNS_FOR_IGSEQ, rp=["I1", "R1", "R2"])

def fastqc_setup_helper_rules():
    """Define per-run rules dynamically to run FastQC at different processing steps"""
    suffix = "_fastqc.quals.csv"
    rps = ["I1", "R1", "R2"]
    for run in RUNS_FOR_IGSEQ:
        samples = [attrs["Sample"] for attrs in SAMPLES.values() if attrs["Run"] == run]
        args={"run": run, "sample": samples, "rp": rps, "suffix": suffix}
        rule:
            f"FastQC for raw R1/R2/I1 reads for run {run}"
            name: f"fastqc_reads_{run}"
            input: expand("analysis/fastqc/reads/{run}/Undetermined_S0_L001_{rp}_001{suffix}", **args)
        rule:
            f"FastQC for demuxed reads for run {run}"
            name: f"fastqc_demux_{run}"
            input: expand("analysis/fastqc/demux/{run}/{sample}.{rp}{suffix}", **args)
        rule:
            f"FastQC for merged reads for run {run}"
            name: f"fastqc_merge_{run}"
            input: expand("analysis/fastqc/merge/{run}/{sample}{suffix}", **args)
fastqc_setup_helper_rules()

rule fastqc:
    output:
        html="analysis/fastqc/{path}_fastqc.html",
        zip="analysis/fastqc/{path}_fastqc.zip"
    input: "analysis/{path}.fastq.gz"
    threads: 8
    shell: "fastqc -t {threads} {input} -o $(dirname {output.html})"

rule fastqc_extract_quals:
    """Extract a FastQC per-base quality summary table as CSV"""
    output: "analysis/fastqc/{path}_fastqc.quals.csv"
    input: "analysis/fastqc/{path}_fastqc.zip"
    shell: "fastqc_extract_quals.py {input} {output}"
