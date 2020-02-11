"""
Custom demultiplexing for our IgSeq protocol.
"""

from igseq.demux import demux

TARGET_DEMUX = expand("demux/{run}/{sample}.{rp}.fastq.gz",
    sample = SAMPLES.keys(),
    run = RUNS.keys(),
    rp = ["R1", "R2", "I1"])

rule all_demux:
    input: TARGET_DEMUX

# This will separate out samples by barcode, and will trim the forward barcode
# from the start of the forward read.  In an ideal world we could use
# cutadapt's features to demultiplex and trim in one go, but it doesn't look
# like it has any option to do this mixed R1/I1 barcoding thing.
rule demux:
    output: expand("demux/{{run}}/{sample}.{rp}.fastq.gz", sample = SAMPLES.keys(), rp = ["R1", "R2", "I1"])
    input:
        r1="data/{run}/Undetermined_S0_L001_R1_001.fastq.gz",
        r2="data/{run}/Undetermined_S0_L001_R2_001.fastq.gz",
        i1="data/{run}/Undetermined_S0_L001_I1_001.fastq.gz"
    params: outdir="demux/{run}"
    log: "logs/demux.{run}.tsv"
    run:
        revcmp = RUNS[wildcards.run]["ReverseComplement"] == "1"
        # Gather samples for just this run
        samples = {}
        for sample_name, sample in SAMPLES.items():
            if sample["Run"] == wildcards.run:
                samples[sample_name] = sample
        with open(log[0], "w") as f_log:
            demux(
                samples=samples,
                fps = {"R1": input.r1, "R2": input.r2, "I1": input.i1},
                outdir=params.outdir, output_files=output,
                revcmp=revcmp, send_stats=f_log)
