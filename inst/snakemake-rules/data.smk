"""
Handlers for data/metadata
"""

from pathlib import Path

PRIMERS = {}
SAMPLES = {}
RUNS = {}
SAMPLES_ALL = set()
try:
    # simple name/seq dict
    PRIMERS = load_primers("metadata/primers.csv")
    # key is run, entries are lists of sample dicts
    SAMPLES = load_samples("metadata/samples.csv", PRIMERS)
    # key is run, entries are dicts
    RUNS = load_runs("metadata/runs.csv")
    # A set of all unique samples
    for samps in SAMPLES.values():
        SAMPLES_ALL = SAMPLES_ALL | set([s["Sample"] for s in samps])
except FileNotFoundError:
    print("Skipping metadata loading; be sure to run get_metadata rule.")
    
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


rule get_metadata:
    output: expand("metadata/{sheet}.csv", sheet = ["samples", "specimens", "primers", "runs"])
    input: "metadata.yml"
    run: R("""
            devtools::load_all("igseq")
            update_metadata_via_yaml("{input}", ".")
            """)
