rule fastqc_merge:
    input: expand("analysis/fastqc/merge/{run}/{sample}_fastqc.html", zip, run=[attrs["Run"] for attrs in SAMPLES.values() if attrs["Run"]], sample=[attrs["Sample"] for attrs in SAMPLES.values() if attrs["Run"]])

rule fastqc_reads:
    input: expand("analysis/fastqc/reads/{run}/Undetermined_S0_L001_{rp}_001_fastqc.html", run=RUNS_FOR_IGSEQ, rp=["I1", "R1", "R2"])

rule fastqc:
    output: "analysis/fastqc/{path}_fastqc.html"
    input: "analysis/{path}.fastq.gz"
    threads: 8
    shell: "fastqc -t {threads} {input} -o $(dirname {output})"
