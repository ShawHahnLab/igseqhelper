"""
Helpers finishing up the per-run processing.
"""

def reads_by_sample_input(wildcards):
    runid = SAMPLES[wildcards.sample]["Run"]
    if not runid:
        raise ValueError("No run given for sample %s" % wildcards.sample)
    return expand(
        "analysis/trim/{run}/{chunk}/{{sample}}.{{rp}}.fastq.gz",
        run=runid, chunk=CHUNKS)

rule reads_by_sample:
    """Aggregate trimmed sequencer samples into per-sample FASTQ files.
    
    This is the end of the road for per-run processing.  We'll drop the
    distinction of what run matches what sample in the paths here so downstream
    rules can just request on a per-sample basis.
    """
    output: "analysis/reads-by-sample/{sample}.{rp}.fastq.gz"
    input: reads_by_sample_input
    run:
        if input:
            shell("cat {input} > {output}")
        else:
            shell("echo -n | gzip > {output}")
