"""
Custom demultiplexing for our IgSeq protocol.
"""

import gzip
from Bio import SeqIO
from igseq.demux import demux, normalize_read_files
from igseq.igblast import tabulate_igblast
from igseq.util import make_chunk_str

def gather_samples_for_demux(samples_all, runid):
    """Make dictionary of sample names to attrs for a given Run ID."""
    samples = {}
    for sample_name, sample in samples_all.items():
        if sample["Run"] == runid:
            samples[sample_name] = sample
    return samples

def demux_input_for(runattrs):
    """Make Snakemake input function for demuxing a given run."""
    keys = ["R1", "R2", "I1"]
    runid=runattrs["Run"]
    if runattrs["DemuxBy"] == "igblast":
        suffix = "fasta"
    elif runattrs["DemuxBy"] == "barcode":
        suffix = "fastq.gz"
    def demux_input(wildcards):
        paths = expand("analysis/{run}/chunk_{chunk}_{rp}.{suffix}",
            run=runid, chunk=wildcards.chunk, rp=keys, suffix=suffix)
        return dict(zip(keys, paths))
    return demux_input

# Inspired by the for loop that makes rules on the fly from IgDiscover:
# https://github.com/NBISweden/IgDiscover/blob/master/src/igdiscover/Snakefile#L387
# Thanks!
# NOTE: be very careful not to reference loop variables directly in the rule
# definition.  Instead save objects as parameters so they're tied to a specific
# iteration.
# This way each rule only needs to create output files specific to that run.
# Before that I had every single sample from every single run, just because of
# how Snakemake works.  Checkpoint rules might be another way to handle this,
# though we *do* know what the outputs will be ahead of time so it wouldn't
# really be the intended use case.
for runid in RUNS.keys():
    # Gather samples for just this run
    samples = gather_samples_for_demux(SAMPLES, runid)
    rule:
        message: "demultiplexing chunk {chunk} for {run}".format(chunk=wildcards.chunk, run=runid)
        output: "analysis/demux/{run}/{chunk}/{sample}.{rp}.fastq.gz", run=runid, sample=samples.keys(), rp=["R1", "R2", "I1"]))
        input: unpack(demux_input_for(runid))
        params:
            samples=samples,
            runattrs=RUNS[runid],
            outdir="analysis/demux/{run}".format(run=runid)
        log: "analysis/logs/demux.{run}.tsv".format(run=runid)
        run:
            with open(log[0], "w") as f_log:
                demux(
                    samples=params.samples,
                    fps=input,
                    runattrs=params.runattrs,
                    outdir=params.outdir,
                    send_stats=f_log)

rule igblast_parse:
    output: "analysis/data/igblast/{run}/chunk_{chunk}.csv"
    input: "analysis/data/igblast/{run}/chunk_{chunk}.txt"
    run: tabulate_igblast(input[0], output[0])

rule igblast_seqids_by_chain:
    """Split AIRR table of blast results in heavy and light."""
    output:
        heavy="analysis/data/igblast/{run}/chunk_{chunk}.heavy.txt",
        light="analysis/data/igblast/{run}/chunk_{chunk}.light.txt"
    input: "analysis/data/igblast/{run}/chunk_{chunk}.txt"
    run: split_igblast_seqids_by_chain(input[0], output.heavy, output.light)

rule igblast_raw:
    output: "analysis/data/igblast/{run}/chunk_{chunk}.txt"
    input:
        query="analysis/data/{run}/chunk_{chunk}_R1.fasta",
        db_v="analysis/data/igblast/v.fasta.nhr",
        db_d="analysis/data/igblast/d.fasta.nhr",
        db_j="analysis/data/igblast/j.fasta.nhr"
    params:
        outfmt=19, # AIRR TSV format
        organism="rhesus_monkey"
    shell:
        """
            dbv={input.db_v}
            dbv=${{dbv%.nhr}}
            dbd={input.db_d}
            dbd=${{dbd%.nhr}}
            dbj={input.db_j}
            dbj=${{dbj%.nhr}}
            igblastn \
                -germline_db_V $dbv \
                -germline_db_D $dbd \
                -germline_db_J $dbj \
                -outfmt {params.outfmt} \
                -organism {params.organism} \
                -ig_seqtype Ig \
                -query {input.query} \
                -out {output}
        """

rule igblast_db:
    output: "analysis/data/igblast/{prefix}.fasta.nhr"
    input: "analysis/data/igblast/{prefix}.fasta"
    shell: "makeblastdb -dbtype nucl -parse_seqids -in {input}"

TARGET_DEMUX = expand(
    outputs_per_run("analysis/demux/{run}/{sample}.{{rp}}.fastq.gz", SAMPLES),
    rp=["R1", "R2", "I1"])

rule all_demux:
    input: TARGET_DEMUX
