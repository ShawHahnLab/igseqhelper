"""
Custom demultiplexing for our IgSeq protocol.
"""

from igseq.demux import demux

# Inspired by the for loop that makes rules on the fly from IgDiscover:
# https://github.com/NBISweden/IgDiscover/blob/master/src/igdiscover/Snakefile#L387
# Thanks!
# NOTE: be very careful not to reference loop variables directly in the rule
# definition.  Instead save objects as parameters so they're tied to a specific
# iteration.
# This way each rule only needs to create output files specific to that run.
# Before that I had every single sample from every single run, just because of
# how Snakemake works.  Checkpoint rules might be another way to handle this.
for runid in RUNS.keys():
    revcmp = RUNS[runid]["ReverseComplement"] == "1"
    # Gather samples for just this run
    samples = {}
    for sample_name, sample in SAMPLES.items():
        if sample["Run"] == runid:
            samples[sample_name] = sample
    # This will separate out samples by barcode, and will trim the forward barcode
    # from the start of the forward read.  In an ideal world we could use
    # cutadapt's features to demultiplex and trim in one go, but it doesn't look
    # like it has any option to do this mixed R1/I1 barcoding thing.
    rule:
        message: "demultiplexing for {run}".format(run=runid)
        output: expand("analysis/demux/{run}/{sample}.{rp}.fastq.gz", run=runid, sample=samples.keys(), rp=["R1", "R2", "I1"])
        input:
            r1="data/{run}/Undetermined_S0_L001_R1_001.fastq.gz".format(run=runid),
            r2="data/{run}/Undetermined_S0_L001_R2_001.fastq.gz".format(run=runid),
            i1="data/{run}/Undetermined_S0_L001_I1_001.fastq.gz".format(run=runid)
        params:
            samples=samples,
            dorevcmp=revcmp,
            outdir="analysis/demux/{run}".format(run=runid)
        log: "analysis/logs/demux.{run}.tsv".format(run=runid)
        run:
            with open(log[0], "w") as f_log:
                demux(
                    samples=params.samples,
                    fps = {"R1": input.r1, "R2": input.r2, "I1": input.i1},
                    outdir=params.outdir, output_files=output,
                    dorevcmp=params.dorevcmp, send_stats=f_log)

TARGET_DEMUX = expand(
    outputs_per_run("analysis/demux/{run}/{sample}.{{rp}}.fastq.gz", SAMPLES),
    rp=["R1", "R2", "I1"])

rule all_demux:
    input: TARGET_DEMUX
