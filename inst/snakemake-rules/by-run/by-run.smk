"""
Rules for dealing with data aggregated at most by a single run.

This includes demultiplexing and adapter trimming.
"""
from igseqhelper.util import make_chunk_str
CHUNKS = make_chunk_str(20)
include: "demux.smk"
include: "trim.smk"
include: "reads.smk"

rule all_reads_by_sample:
    input: expand("analysis/reads-by-sample/{sample}.{rp}.fastq.gz", sample=SAMPLES.keys(), rp=["R1", "R2"])
