"""
pRESTO rules for our protocol.

I needed to use our own custom demultiplexing first and change the defaults in
a few places to account for our weird method.

See also:

 * https://github.com/ressy/example-presto
 * https://presto.readthedocs.io/en/stable/workflows/Greiff2014_Workflow.html
"""

from igseq.presto import (PRESTO_OPTS, prep_primers_fwd)

TARGET_PRESTO_DATA = expand("presto/data/{sample}_1.fastq", sample=SAMPLES.keys())
TARGET_PRESTO_ASSEMBLY = expand("presto/assemble/{sample}_assemble-pass.fastq", sample=SAMPLES.keys())
TARGET_PRESTO_QC = expand("presto/qual/{sample}-FWD_primers-pass.fastq", sample=SAMPLES.keys())
TARGET_PRESTO_ALL = expand("presto/collapse/{sample}_atleast-2.fastq", sample=SAMPLES.keys())

rule all_presto_data:
    input: TARGET_PRESTO_DATA

rule all_presto_assembly:
    input: TARGET_PRESTO_ASSEMBLY

rule all_presto_qc:
    input: TARGET_PRESTO_QC

rule all_presto:
    input: TARGET_PRESTO_ALL

### Before pRESTO: demultiplex, trim reads, and gather input data

rule presto_data:
    output:
        r1="presto/data/{sample}_1.fastq",
        r2="presto/data/{sample}_2.fastq"
    input:
        r1=expand("trim/{run}/{{sample}}.R1.fastq.gz", run = RUNS.keys()),
        r2=expand("trim/{run}/{{sample}}.R2.fastq.gz", run = RUNS.keys())
    shell:
        """
            zcat {input.r1} > {output.r1}
            zcat {input.r2} > {output.r2}
        """

rule presto_primers:
    output: fwd="presto/vprimers.fasta"
    input: sequences="metadata/sequences.csv"
    run: prep_primers_fwd(input.sequences, output.fwd)

### Paired-end Assembly

rule presto_assembly:
    output: "presto/assemble/{sample}_assemble-pass.fastq"
    input:
        r1="presto/data/{sample}_1.fastq",
        r2="presto/data/{sample}_2.fastq"
    log: "logs/presto/assemble.{sample}.log"
    params:
        coord=PRESTO_OPTS["assembly"]["coord"],
        rc=PRESTO_OPTS["assembly"]["rc"],
        minlen=PRESTO_OPTS["assembly"]["minlen"],
        maxlen=PRESTO_OPTS["assembly"]["maxlen"]
    threads: 8
    shell:
        """
            AssemblePairs.py align \
                --minlen {params.minlen} --maxlen {params.maxlen} \
                -1 {input.r1} -2 {input.r2} \
                --nproc {threads} -o {output} --log {log} \
                --coord {params.coord} --rc {params.rc}
        """

### Quality Control

rule presto_qual:
    output: "presto/qual/{sample}_quality-pass.fastq"
    input: "presto/assemble/{sample}_assemble-pass.fastq"
    log: "logs/presto/qual.{sample}.log"
    params:
        mean_qual=PRESTO_OPTS["qc"]["mean_qual"]
    threads: 8
    shell:
        """
            FilterSeq.py quality \
                --nproc {threads} -s {input} -o {output} --log {log} \
                -q {params.mean_qual}
                
        """

rule presto_qual_fwd:
    output: "presto/qual/{sample}-FWD_primers-pass.fastq"
    input:
        data="presto/qual/{sample}_quality-pass.fastq",
        primers="presto/vprimers.fasta"
    log: "logs/presto/qual_fwd.{sample}.log"
    params:
        start=PRESTO_OPTS["qc"]["fwd_start"],
        mode=PRESTO_OPTS["qc"]["fwd_mode"],
        pf=PRESTO_OPTS["qc"]["fwd_pf"]
    threads: 8
    shell:
        """
            MaskPrimers.py score \
                --nproc {threads} -s {input.data} -o {output} --log {log} \
                -p {input.primers} \
                --start {params.start} --mode {params.mode} --pf {params.pf}
        """

# Skipping this.
#
# I don't think we have any trace of the reverse primer visible in the sequence
# itself, only in the I1 read.
rule presto_qual_rev:
    output: "presto/qual/{sample}-REV_primers-pass.fastq"
    input:
        data="presto/qual/{sample}-FWD_primers-pass.fastq",
        primers=""
    log: "logs/presto/qual_rev.{sample}.log"
    params:
        start=PRESTO_OPTS["qc"]["rev_start"],
        mode=PRESTO_OPTS["qc"]["rev_mode"],
        pf=PRESTO_OPTS["qc"]["rev_pf"]
    threads: 8
    shell:
        """
            MaskPrimers.py score \
                    --nproc {threads} -s {input.data} -o {output} --log {log} \
                    -p {input.primers} \
                    --start {params.start} --mode {params.mode} --pf {params.pf} --revpr
        """

### Deduplication and Filtering

rule presto_collapse:
    output: "presto/collapse/{sample}_collapse-unique.fastq"
    input: "presto/qual/{sample}-FWD_primers-pass.fastq"
    log: "logs/presto/collapse.{sample}.log"
    params:
        uf=PRESTO_OPTS["collapse"]["uf"],
        cf=PRESTO_OPTS["collapse"]["cf"]
    shell:
        """
            CollapseSeq.py \
                -s {input} -o {output} --log {log} \
                -n 20 --inner --uf {params.uf} \
                --cf {params.cf} --act set \
        """

rule presto_collapse_atleast:
    output: "presto/collapse/{sample}_atleast-2.fastq"
    input: "presto/collapse/{sample}_collapse-unique.fastq"
    params:
        f=PRESTO_OPTS["collapse"]["f"],
        num=PRESTO_OPTS["collapse"]["num"]
    shell:
        """
            SplitSeq.py group \
                -s {input} --outname {wildcards.sample} \
                -f {params.f} --num {params.num}
        """
