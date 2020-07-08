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

from igseq.data import amplicon_files
from igseq.presto import (PRESTO_OPTS, prep_primers_fwd)

TARGET_PRESTO_DATA = expand(amplicon_files("analysis/presto/data/{chain}.{chain_type}/{specimen}.{{rp}}.fastq", SAMPLES, "IgG+"), rp=["R1", "R2"])
TARGET_PRESTO_ASSEMBLY = amplicon_files("analysis/presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq", SAMPLES, "IgG+")
TARGET_PRESTO_QC = amplicon_files("analysis/presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq", SAMPLES, "IgG+")
TARGET_PRESTO_ALL = amplicon_files("analysis/presto/collapse/{chain}.{chain_type}/{specimen}_atleast-2.fastq", SAMPLES, "IgG+")

rule all_presto_data:
    input: TARGET_PRESTO_DATA

rule all_presto_assembly:
    input: TARGET_PRESTO_ASSEMBLY

rule all_presto_qc:
    input: TARGET_PRESTO_QC

rule all_presto:
    input: TARGET_PRESTO_ALL

### Before pRESTO: demultiplex, trim reads. now, gather input data

rule presto_data:
    """Aggregate trimmed sequencer samples into per-specimen FASTQ files.

    pRESTO uses plaintext FASTQ so we'll leave them uncompressed, too.
    """
    output: "analysis/presto/data/{chain}.{chain_type}/{specimen}.{rp}.fastq"
    input: "analysis/reads-by-specimen/{chain}.{chain_type}/{specimen}.{rp}.fastq.gz"
    shell: "zcat {input} > {output}"

rule presto_primers:
    output: fwd="analysis/presto/vprimers.fasta"
    input: sequences="metadata/sequences.csv"
    run: prep_primers_fwd(input.sequences, output.fwd)

### Paired-end Assembly

rule presto_assembly:
    output: "analysis/presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq"
    input:
        r1="analysis/presto/data/{chain}.{chain_type}/{specimen}.R1.fastq",
        r2="analysis/presto/data/{chain}.{chain_type}/{specimen}.R2.fastq"
    log: "analysis/logs/presto/assemble.{chain}.{chain_type}.{specimen}.log"
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
    output: "analysis/presto/qual/{chain}.{chain_type}/{specimen}_quality-pass.fastq"
    input: "analysis/presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq"
    log: "analysis/logs/presto/qual.{chain}.{chain_type}.{specimen}.log"
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
    output: "analysis/presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq"
    input:
        data="analysis/presto/qual/{chain}.{chain_type}/{specimen}_quality-pass.fastq",
        primers="analysis/presto/vprimers.fasta"
    log: "analysis/logs/presto/qual_fwd.{chain}.{chain_type}.{specimen}.log"
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
    output: "analysis/presto/qual/{chain}.{chain_type}/{specimen}-REV_primers-pass.fastq"
    input:
        data="analysis/presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq",
        primers=""
    log: "analysis/logs/presto/qual_rev.{chain}.{chain_type}.{specimen}.log"
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
    output: "analysis/presto/collapse/{chain}.{chain_type}/{specimen}_collapse-unique.fastq"
    input: "analysis/presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq"
    log: "analysis/logs/presto/collapse.{chain}.{chain_type}.{specimen}.log"
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
    output: "analysis/presto/collapse/{chain}.{chain_type}/{specimen}_atleast-2.fastq"
    input: "analysis/presto/collapse/{chain}.{chain_type}/{specimen}_collapse-unique.fastq"
    params:
        f=PRESTO_OPTS["collapse"]["f"],
        num=PRESTO_OPTS["collapse"]["num"]
    shell:
        """
            SplitSeq.py group \
                -s {input} --outname {wildcards.specimen} \
                -f {params.f} --num {params.num}
        """
