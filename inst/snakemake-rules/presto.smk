"""
pRESTO rules for our protocol.

I needed to use our own custom demultiplexing first and change the defaults in
a few places to account for our weird method.  At this point individual
sequencer samples are aggregated together on a specimen and cell type basis, so
that duplicate sequences for the same phsyical specimen are correctly handled
as such.

See also:

 * https://github.com/ressy/example-presto
 * https://presto.readthedocs.io/en/stable/workflows/Greiff2014_Workflow.html
"""

from igseq.presto import (PRESTO_OPTS, prep_primers_fwd, specimens_per_sample)

TARGET_PRESTO_DATA = expand(specimens_per_sample("presto/data/{chain}.{chain_type}/{specimen}.{{rp}}.fastq", SAMPLES), rp=["R1", "R2"])
TARGET_PRESTO_ASSEMBLY = expand(specimens_per_sample("presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq", SAMPLES))
TARGET_PRESTO_QC = expand(specimens_per_sample("presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq", SAMPLES))
TARGET_PRESTO_ALL = expand(specimens_per_sample("presto/collapse/{chain}.{chain_type}/{specimen}_atleast-2.fastq", SAMPLES))

rule all_presto_data:
    input: TARGET_PRESTO_DATA

rule all_presto_assembly:
    input: TARGET_PRESTO_ASSEMBLY

rule all_presto_qc:
    input: TARGET_PRESTO_QC

rule all_presto:
    input: TARGET_PRESTO_ALL

### Before pRESTO: demultiplex, trim reads. now, gather input data

def trimmed_specimen_samples(wildcards):
    """Get trimmed sequencer samples matching specimen+antibody type.

    We may have prepped and sequenced the same specimen multiple times so this
    will give a list of one or more files.
    """
    pattern = "trim/{run}/{sample}.{rp}.fastq.gz"
    filepaths = []
    for samp_name, samp_attrs in SAMPLES.items():
        if \
            samp_attrs["Specimen"] == wildcards.specimen and \
            samp_attrs["Chain"] == wildcards.chain and \
            samp_attrs["Type"] == wildcards.chain_type:
            filepaths.append(pattern.format(
                run=samp_attrs["Run"],
                sample=samp_name,
                rp=wildcards.rp))
    return filepaths

rule presto_data:
    """Aggregate trimmed sequencer samples into per-specimen FASTQ files.

    pRESTO uses plaintext FASTQ so we'll leave them uncompressed, too.
    """
    output: "presto/data/{chain}.{chain_type}/{specimen}.{rp}.fastq"
    input: trimmed_specimen_samples
    shell: "zcat {input} > {output}"

rule presto_primers:
    output: fwd="presto/vprimers.fasta"
    input: sequences="metadata/sequences.csv"
    run: prep_primers_fwd(input.sequences, output.fwd)

### Paired-end Assembly

rule presto_assembly:
    output: "presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq"
    input:
        r1="presto/data/{chain}.{chain_type}/{specimen}.R1.fastq",
        r2="presto/data/{chain}.{chain_type}/{specimen}.R2.fastq"
    log: "logs/presto/assemble.{chain}.{chain_type}.{specimen}.log"
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
    output: "presto/qual/{chain}.{chain_type}/{specimen}_quality-pass.fastq"
    input: "presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq"
    log: "logs/presto/qual.{chain}.{chain_type}.{specimen}.log"
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
    output: "presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq"
    input:
        data="presto/qual/{chain}.{chain_type}/{specimen}_quality-pass.fastq",
        primers="presto/vprimers.fasta"
    log: "logs/presto/qual_fwd.{chain}.{chain_type}.{specimen}.log"
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
    output: "presto/qual/{chain}.{chain_type}/{specimen}-REV_primers-pass.fastq"
    input:
        data="presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq",
        primers=""
    log: "logs/presto/qual_rev.{chain}.{chain_type}.{specimen}.log"
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
    output: "presto/collapse/{chain}.{chain_type}/{specimen}_collapse-unique.fastq"
    input: "presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq"
    log: "logs/presto/collapse.{chain}.{chain_type}.{specimen}.log"
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
    output: "presto/collapse/{chain}.{chain_type}/{specimen}_atleast-2.fastq"
    input: "presto/collapse/{chain}.{chain_type}/{specimen}_collapse-unique.fastq"
    params:
        f=PRESTO_OPTS["collapse"]["f"],
        num=PRESTO_OPTS["collapse"]["num"]
    shell:
        """
            SplitSeq.py group \
                -s {input} --outname {wildcards.specimen} \
                -f {params.f} --num {params.num}
        """
