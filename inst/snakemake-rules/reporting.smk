### Reporting

rule fastq_counts:
    output: "{prefix}.fastq.gz.counts"
    input: "{prefix}.fastq.gz"
    shell: "zcat {input} | wc -l | awk -v OFS=, -v f={input} '{{print f,$0/4}}' > {output}"

rule fasta_counts:
    output: "{prefix}.fasta.counts"
    input: "{prefix}.fasta"
    shell: "grep -c '^>' {input} | awk -v OFS=, -v f={input} '{{print f,$0}}' > {output}"

rule counts_table:
    output: "counts.csv"
    input:
        raw=expand("0-data/{run}/" + RAW + ".counts", run = RUNS.keys(), rp = ["R1"])
    shell:
        """
            echo "Filename,NumSequences" > {output}
            cat {input} | sort >> {output}
        """
