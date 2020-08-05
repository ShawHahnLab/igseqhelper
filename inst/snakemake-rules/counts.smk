"""
Generic helper rules to get sequence counts from various formats.

These use Snakemake's support for path separators within wildcards to match
arbitrary file paths, and output to a separate directory tree.
"""

rule fqgz_read_counts:
    output: "analysis/counts/{path}.fastq.gz.counts"
    input: "analysis/{path}.fastq.gz"
    shell: "zcat {input} | sed -n 1~4p | wc -l > {output}"

rule fastq_read_counts:
    output: "analysis/counts/{path}.fastq.counts"
    input: "analysis/{path}.fastq"
    shell: "sed -n 1~4p {input} | wc -l > {output}"

rule fagz_read_counts:
    output: "analysis/counts/{path}.fasta.gz.counts"
    input: "analysis/{path}.fasta.gz"
    shell: "zgrep -c '^>' {input} > {output}"

rule fasta_read_counts:
    output: "analysis/counts/{path}.fasta.counts"
    input: "analysis/{path}.fasta"
    shell: "grep -c '^>' {input} > {output}"
