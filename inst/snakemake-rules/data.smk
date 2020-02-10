"""
Handlers for data/metadata
"""

from pathlib import Path

def _setup_metadata(fp_primers, fp_samples, fp_runs):
    global SEQUENCES, SAMPLES, RUNS, SAMPLES_ALL
    # key is sequence name, entries are dicts
    SEQUENCES = load_sequences("metadata/sequences.csv")
    # key is run, entries are lists of sample dicts
    SAMPLES = load_samples("metadata/samples.csv", SEQUENCES)
    # key is run, entries are dicts
    RUNS = load_runs("metadata/runs.csv")
    # A set of all unique samples
    SAMPLES_ALL = set()
    for samps in SAMPLES.values():
        SAMPLES_ALL = SAMPLES_ALL | set([s["Sample"] for s in samps])

try:
    _setup_metadata("metadata/sequences.csv", "metadata/samples.csv", "metadata/runs.csv")
except FileNotFoundError:
    print("Skipping metadata loading; be sure to run get_metadata rule.")
    SEQUENCES = {}
    SAMPLES = {}
    RUNS = {}
    SAMPLES_ALL = set()

SAMPLES_H = []
SAMPLES_L = []
ALLELES = {"heavy": {"V": "", "D": "", "J": ""}, "light": {"V": "", "J": ""}}

RAW = "Undetermined_S0_L001_{rp}_001.fastq.gz"

def get_data(runid, outdir):
    """Get the raw data files for a single run, either from local disk or from URLs."""
    LOGGER.info("get_data: runid %s", runid)
    LOGGER.info("get_data: outdir %s", outdir)
    rundir = Path("/seq/runs" / runid)
    if rundir.is_dir():
        LOGGER.info("get_data: running bcl2fastq")
        shell(
            """
                bcl2fastq --create-fastq-for-index-reads \
                    -R /seq/runs/%s \
                    -o %s
            """ % (runid, outdir))
    else:
        url = RUNS[runid].get("URL")
        if url:
            LOGGER.info("get_data: downloading single zip")
            shell("cd data/{wildcards.run} && wget '%s' && unzip {wildcards.run}.zip" % url)
        else:
            url_r1 = RUNS[runid].get("URLR1")
            url_r2 = RUNS[runid].get("URLR2")
            url_i1 = RUNS[runid].get("URLI1")
            if url_r1 and url_r2 and url_r3:
                LOGGER.info("get_data: downloading separate zips")
                shell(
                    """
                        cd data/{wildcards.run}
                        wget '%s' && unzip %s_R1.zip
                        wget '%s' && unzip %s_R2.zip
                        wget '%s' && unzip %s_I1.zip
                    """ % (url_r1, runid, url_r2, runid, url_i1, runid))
            else:
                raise ValueError("Need raw data or URLs for run %s" % runid)

rule all_get_data:
    input:
        files=expand("data/{run}/" + RAW, run = SAMPLES.keys(), rp = ["R1", "R2", "I1"]),
        metadata=lambda w: checkpoints.get_metadata.get().output

rule get_data:
    """Get data files for a run.

    If raw data is available locally, run Illumina's bcl2fastq to create R1/R2/I1
    files.  (Add location of bcl2fatq to PATH if needed.)  If not, download from
    URLs listed in metadata.
    """
    output: expand("data/{run}/" + RAW, run = "{run}", rp = ["R1", "R2", "I1"])
    input: lambda w: checkpoints.get_metadata.get().output
    run: get_data(wildcards.run, Path(output[0]).parent)

checkpoint get_metadata:
    output:
       samples="metadata/samples.csv",
       specimens="metadata/specimens.csv",
       sequences="metadata/sequences.csv",
       runs="metadata/runs.csv"
    input: "metadata.yml"
    run:
        R("""
          devtools::load_all("igseq")
          update_metadata_via_yaml("{input}", ".")
          """)
        _setup_metadata(input.sequences, input.samples, input.runs)
