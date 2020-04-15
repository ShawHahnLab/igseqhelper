"""
Sample-specific adapter trimming.

Our barcodes vary in length, so we'll call cutadapt with sample-specific
barcode arguments.
"""

from igseq.trim import (adapter_fwd, adapter_rev)

TARGET_TRIMMED = expand(
    outputs_per_run("analysis/trim/{run}/{{chunk}}/{sample}.{{rp}}.fastq.gz", SAMPLES),
    chunk=CHUNKS, rp=["R1", "R2"])

rule all_trim:
    input: TARGET_TRIMMED

# This will trim off everything from the ends of R1 and R2 that isn't part of
# the biological sequence (so, barcodes + primers) and will apply a slight
# quality threshold.  This ensures that pRESTO's alignment step has the best
# chance possible of properly pairing reads and that no technical bits have any
# chance of getting mixed in with the biological sequence content.
rule trim:
    output:
        r1="analysis/trim/{run}/{chunk}/{sample}.R1.fastq.gz",
        r2="analysis/trim/{run}/{chunk}/{sample}.R2.fastq.gz"
    input:
        r1="analysis/demux/{run}/{chunk}/{sample}.R1.fastq.gz",
        r2="analysis/demux/{run}/{chunk}/{sample}.R2.fastq.gz"
    threads: 4
    params:
        adapter_R1=lambda w: adapter_fwd(SAMPLES[w.sample], SEQUENCES),
        adapter_R2=lambda w: adapter_rev(SAMPLES[w.sample], SEQUENCES),
        # https://cutadapt.readthedocs.io/en/stable/guide.html#filtering-reads
        # PEAR may have trouble with zero-length reads, I think
        min_len=50,
        # https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming
        # https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm
        quality_cutoff=15
    log: "analysis/logs/trim.{run}.{chunk}.{sample}.log"
    shell:
        """
            cutadapt --cores {threads} --quality-cutoff {params.quality_cutoff} \
                -m {params.min_len} \
                -a {params.adapter_R1} -A {params.adapter_R2} \
                -o {output.r1} -p {output.r2} \
                {input.r1} {input.r2} > {log}
        """
