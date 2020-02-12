"""
Generic helper rules to get sequence counts from various formats.

These use Snakemake's support for path separators within wildcards to match
arbitrary file paths, and output to a separate directory tree.
"""

rule fqgz_read_counts:
    output: "counts/{path}.fastq.gz.counts"
    input: "{path}.fastq.gz"
    shell: "zcat {input} | sed -n 1~4p | wc -l > {output}"

rule fastq_read_counts:
    output: "counts/{path}.fastq.counts"
    input: "{path}.fastq"
    shell: "sed -n 1~4p {input} | wc -l > {output}"

rule fagz_read_counts:
    output: "counts/{path}.fasta.gz.counts"
    input: "{path}.fasta.gz"
    shell: "zgrep -c '^>' {input} > {output}"

rule fasta_read_counts:
    output: "counts/{path}.fasta.counts"
    input: "{path}.fasta"
    shell: "grep -c '^>' {input} > {output}"
