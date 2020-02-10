"""
Handlers for data/metadata
"""

from pathlib import Path

def _setup_metadata(fp_primers, fp_samples, fp_runs):
    global PRIMERS SAMPLES RUNS SAMPLES_ALL
    # simple name/seq dict
    PRIMERS = load_primers("metadata/primers.csv")
    # key is run, entries are lists of sample dicts
    SAMPLES = load_samples("metadata/samples.csv", PRIMERS)
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
    PRIMERS = {}
    SAMPLES = {}
    RUNS = {}
    SAMPLES_ALL = set()

SAMPLES_H = []
SAMPLES_L = []
ALLELES = {"heavy": {"V": "", "D": "", "J": ""}, "light": {"V": "", "J": ""}}

RAW = "Undetermined_S0_L001_{rp}_001.fastq.gz"

rule all_get_data:
    input: expand("data/{run}/" + RAW, run = SAMPLES.keys(), rp = ["R1", "R2", "I1"])

# If raw data is available locally, run Illumina's bcl2fastq to create R1/R2/I1
# files.  (Add location of bcl2fatq to PATH if needed.)  If not, download from
# URLs listed in metadata.
rule get_data:
    output: expand("data/{run}/" + RAW, run = "{run}", rp = ["R1", "R2", "I1"])
    input: lambda w: checkpoints.get_metadata.get().output
    run:
        rundir = Path("/seq/runs/{wildcards.run}")
        if rundir.is_dir():
            shell(
                """
                    bcl2fastq --create-fastq-for-index-reads \
                        -R /seq/runs/{wildcards.run} \
                        -o $(dirname {output[0]})
                """)
        else:
            url = RUNS[wildcards.run].get("URL")
            if url:
                shell("cd data/{wildcards.run} && wget '%s' && unzip {wildcards.run}.zip" % url)
            else:
                url_r1 = RUNS[wildcards.run].get("URLR1")
                url_r2 = RUNS[wildcards.run].get("URLR2")
                url_i1 = RUNS[wildcards.run].get("URLI1")
                if url_r1 and url_r2 and url_r3:
                    shell(
                        """
                            cd data/{wildcards.run}
                            wget '%s' && unzip {wildcards.run}_R1.zip
                            wget '%s' && unzip {wildcards.run}_R2.zip
                            wget '%s' && unzip {wildcards.run}_I1.zip
                        """ % (url_r1, url_r2, url_i1))
                else: 
                    raise ValueError("Need raw data or URLs for run %s" % wildcards.run)


checkpoint get_metadata:
    output:
       samples="metadata/samples.csv"
       specimens="metadata/specimens.csv"
       sequences="metadata/sequences.csv"
       runs="metadata/runs.csv"
    input: "metadata.yml"
    run:
        R("""
          devtools::load_all("igseq")
          update_metadata_via_yaml("{input}", ".")
          """)
        _setup_metadata(input.sequences, input.samples, input.runs)
