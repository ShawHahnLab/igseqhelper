"""
Sample-specific adapter trimming.

Our barcodes vary in length, so we'll call cutadapt with sample-specific
barcode arguments.
"""

from igseq.trim import (adapter_fwd, adapter_rev)

TARGET_TRIMMED = expand(
    outputs_per_run("trim/{run}/{sample}.{{rp}}.fastq.gz", SAMPLES),
    rp=["R1", "R2"])

rule all_trim:
    input: TARGET_TRIMMED

# This will trim off everything from the ends of R1 and R2 that isn't part of
# the biological sequence (so, barcodes + primers) and will apply a slight
# quality threshold.  This ensures that pRESTO's alignment step has the best
# chance possible of properly pairing reads and that no technical bits have any
# chance of getting mixed in with the biological sequence content.
rule trim:
    output:
        r1="trim/{run}/{sample}.R1.fastq.gz",
        r2="trim/{run}/{sample}.R2.fastq.gz"
    input:
        r1="demux/{run}/{sample}.R1.fastq.gz",
        r2="demux/{run}/{sample}.R2.fastq.gz"
    threads: 4
    params:
        adapter_R1=lambda w: adapter_fwd(SAMPLES[w.sample], SEQUENCES),
        adapter_R2=lambda w: adapter_rev(SAMPLES[w.sample], SEQUENCES),
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
