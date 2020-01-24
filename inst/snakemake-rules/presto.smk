"""
pRESTO rules for our protocol.

I needed to use our own custom demultiplexing first and change the defaults in
a few places to account for our weird method.

See also:

 * https://github.com/ressy/example-presto
 * https://presto.readthedocs.io/en/stable/workflows/Greiff2014_Workflow.html
"""

from igseq.demux import demux
from igseq.trim import (adapter_fwd, adapter_rev)
from igseq.presto import (PRESTO_OPTS, prep_primers_fwd)

TARGET_PRESTO_DEMUX = expand("presto/demux/{run}/{sample}.{rp}.fastq.gz",
    sample = SAMPLES_ALL,
    run = RUNS.keys(),
    rp = ["R1", "R2", "I1"])
TARGET_PRESTO_TRIMMED = expand("presto/trim/{run}/{sample}.{rp}.fastq.gz",
    sample = SAMPLES_ALL,
    run = RUNS.keys(),
    rp = ["R1", "R2"])
TARGET_PRESTO_DATA = expand("presto/data/{sample}_1.fastq", sample=SAMPLES_ALL)
TARGET_PRESTO_ASSEMBLY = expand("presto/assemble/{sample}_assemble-pass.fastq", sample=SAMPLES_ALL)
TARGET_PRESTO_QC = expand("presto/qual/{sample}-FWD_primers-pass.fastq", sample=SAMPLES_ALL)
TARGET_PRESTO_ALL = expand("presto/collapse/{sample}_atleast-2.fastq", sample=SAMPLES_ALL)

rule all_presto_demux:
    input: TARGET_PRESTO_DEMUX

rule all_presto_trim:
    input: TARGET_PRESTO_TRIMMED

rule all_presto_data:
    input: TARGET_PRESTO_DATA

rule all_presto_assembly:
    input: TARGET_PRESTO_ASSEMBLY

rule all_presto_qc:
    input: TARGET_PRESTO_QC

rule all_presto:
    input: TARGET_PRESTO_ALL

### Before pRESTO: demultiplex, trim reads, and gather input data

# This will separate out samples by barcode, and will trim the forward barcode
# from the start of the forward read.  In an ideal world we could use
# cutadapt's features to demultiplex and trim in one go, but it doesn't look
# like it has any option to do this mixed R1/I1 barcoding thing.
rule presto_demux:
    output: expand("presto/demux/{{run}}/{sample}.{rp}.fastq.gz", sample = SAMPLES_ALL, rp = ["R1", "R2", "I1"])
    input:
        r1="data/{run}/Undetermined_S0_L001_R1_001.fastq.gz",
        r2="data/{run}/Undetermined_S0_L001_R2_001.fastq.gz",
        i1="data/{run}/Undetermined_S0_L001_I1_001.fastq.gz"
    params: outdir="presto/demux/{run}"
    log: "logs/demux.{run}.tsv"
    run:
        revcmp = RUNS[wildcards.run]["ReverseComplement"] == "1"
        with open(log[0], "w") as f_log:
            demux(
                samples = SAMPLES[wildcards.run],
                fps = {"R1": input.r1, "R2": input.r2, "I1": input.i1},
                outdir=params.outdir, output_files=output,
                revcmp=revcmp, send_stats=f_log)

# This will trim off everything from the ends of R1 and R2 that isn't part of
# the biological sequence (so, barcodes + primers) and will apply a slight
# quality threshold.  This ensures that pRESTO's alignment step has the best
# chance possible of properly pairing reads and that no technical bits have any
# chance of getting mixed in with the biological sequence content.
rule presto_trim:
    output:
        r1="presto/trim/{run}/{sample}.R1.fastq.gz",
        r2="presto/trim/{run}/{sample}.R2.fastq.gz"
    input:
        r1="presto/demux/{run}/{sample}.R1.fastq.gz",
        r2="presto/demux/{run}/{sample}.R2.fastq.gz"
    threads: 4
    params:
        adapter_R1=lambda w: adapter_fwd(w.sample, SAMPLES),
        adapter_R2=lambda w: adapter_rev(w.sample, SAMPLES),
        # https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming
        # https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm
        quality_cutoff=10
    log: "logs/trim.{run}.{sample}.log"
    shell:
        """
            cutadapt --cores {threads} --quality-cutoff {params.quality_cutoff} \
                -a {params.adapter_R1} -A {params.adapter_R2} \
                -o {output.r1} -p {output.r2} \
                {input.r1} {input.r2} > {log}
        """

rule presto_data:
    output:
        r1="presto/data/{sample}_1.fastq",
        r2="presto/data/{sample}_2.fastq"
    input:
        r1=expand("presto/trim/{run}/{{sample}}.R1.fastq.gz", run = RUNS.keys()),
        r2=expand("presto/trim/{run}/{{sample}}.R2.fastq.gz", run = RUNS.keys())
    shell:
        """
            zcat {input.r1} > {output.r1}
            zcat {input.r2} > {output.r2}
        """

rule presto_primers:
    output: fwd="presto/vprimers.fasta"
    input: primers="metadata/primers.csv"
    run: prep_primers_fwd(input.primers, output.fwd)

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
