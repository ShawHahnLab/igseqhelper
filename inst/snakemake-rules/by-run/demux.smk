"""
Custom demultiplexing for our IgSeq protocol.
"""

from pathlib import Path
from igseq.demux import demux
from igseq.util import normalize_read_files

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
    suffix = "fastq.gz"
    def demux_input(wildcards):
        paths = expand("analysis/data/{run}/chunk_{chunk}_{rp}.{suffix}",
            run=runid, chunk=wildcards.chunk, rp=keys, suffix=suffix)
        targets = dict(zip(keys, paths))
        return targets
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
        message: "demultiplexing chunk {{wildcards.chunk}} for {run}".format(run=runid)
        output: expand("analysis/demux/{run}/{{chunk}}/{sample}.{rp}.fastq.gz", run=runid, sample=samples.keys(), rp=["R1", "R2", "I1"])
        input: unpack(demux_input_for(RUNS[runid]))
        params:
            samples=samples,
            runattrs=RUNS[runid],
        log: "analysis/logs/demux/{run}/{{chunk}}.tsv".format(run=runid)
        run:
            with open(log[0], "w") as f_log:
                demux(
                    samples=params.samples,
                    fps=dict(input),
                    runattrs=params.runattrs,
                    outdir=Path(output[0]).parent,
                    send_stats=f_log)

rule chunk_raw_data:
    """Split all fastq.gz files from run into fixed number of chunks.

    This will let us parallelize the demultiplexing and any intermediate steps
    and hides inconsistencies between raw outputs between runs (some maybe were
    purportedly demultiplexed by the sequencer, some not).
    """
    output: expand("analysis/data/{{run}}/chunk_{chunk}_{{rp}}.fastq.gz", chunk=CHUNKS)
    input: "data/{run}"
    run:
        fp_sets = normalize_read_files(input[0])
        input_paths = [x[wildcards.rp] for x in fp_sets]
        igseq.data.chunk_fqgz(input_paths, output)

rule igblast_raw:
    output: "analysis/data/igblast/{run}/chunk_{chunk}_airr.tsv"
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
    outputs_per_run("analysis/demux/{run}/{{chunk}}/{sample}.{{rp}}.fastq.gz", SAMPLES),
    chunk=CHUNKS,
    rp=["R1", "R2", "I1"])

rule all_demux:
    input: TARGET_DEMUX
