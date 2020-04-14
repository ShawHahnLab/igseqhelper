"""
Handlers for data.

Raw data initially shows up in data/{run}/ and then has early preprocessing
done in analysis/data/{run}/.
"""

from pathlib import Path
import igseq.data
from igseq.util import make_chunk_str, normalize_read_files

rule all_get_data:
    """Get data files for all sequencing runs.

    This lists the metadata checkpoint output as input so that we know what the
    sequencing runs are before running the rule.
    """
    input:
        rundir=expand("data/{run}", run = SAMPLES.keys()),
        metadata=lambda w: checkpoints.get_metadata.get().output

rule raw_fasta:
    """Convert raw fastq.gz into FASTA for tools that insist on that."""
    output: "analysis/data/{run}/chunk_{chunk}_{rp}.fasta"
    input: "analysis/data/{run}/chunk_{chunk}_{rp}.fastq.gz"
    shell: "seqmagick convert {input} {output}"

rule chunk_raw_data:
    """Split all fastq.gz files from run into fixed number of chunks.

    This will let us parallelize the demultiplexing and any intermediate steps
    and hides inconsistencies between raw outputs between runs (some maybe were
    purportedly demultiplexed by the sequencer, some not).
    """
    output: expand("analysis/data/{{run}}/chunk_{chunk}_{{rp}}.fastq.gz", chunk=make_chunk_str(20))
    input: "data/{run}"
    run:
        fp_sets = normalize_read_files(input[0])
        input_paths = [x[wildcards.rp] for x in fp_sets]
        igseq.data.chunk_fqgz(input_paths, output)

rule get_data:
    """Get data files for a run.

    If raw data is available locally, run Illumina's bcl2fastq to create R1/R2/I1
    files.  (Add location of bcl2fatq to PATH if needed.)  If not, download from
    URLs listed in metadata.  The files in the data/{run}/ directory should then
    be left as-is with everything else happening under analysis/.

    The output here is the run directory itself, since there are a varying
    number of ouputs that can be used in various ways.  An alternative would be
    to use checkpoints to be able to run on a per-file basis, but we generally
    have no use for the distinctions between the separate files (except for
    what's R1, R2, and I1).
    """
    output: protected(directory("data/{run}"))
    input: lambda w: checkpoints.get_metadata.get().output
    params: runpath="/seq/runs"
    run: igseq.data.get_data(wildcards.run, Path(output[0]).parent, RUNS, params.runpath)

#############

def outputs_per_run(pattern, samples):
    """Helper to produce output filenames specific to run/sample combos.

    Give a string pattern with {run} and {sample} and the dictionary of sample
    information and the right run/sample combinations will be filled in.  Be
    sure to escape any other template variables like {{rp}}.
    """
    samples_per_run = igseq.data.get_samples_per_run(samples)
    target = []
    for runid in samples_per_run.keys():
        target.extend(expand(
            pattern,
            sample=samples_per_run[runid],
            run=runid))
    return target
