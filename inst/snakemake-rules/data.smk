"""
Handlers for data/metadata
"""

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

# Run Illumina's bcl2fastq to create R1/R2/I1 files.  Add location of bcl2fatq
# to PATH if needed.
rule get_data:
    output: expand("data/{run}/" + RAW, run = "{run}", rp = ["R1", "R2", "I1"])
    shell:
        """
            bcl2fastq --create-fastq-for-index-reads \
                -R /seq/runs/{wildcards.run} \
                -o $(dirname {output[0]})
        """

rule get_metadata:
    output: expand("metadata/{sheet}.csv", sheet = ["samples", "specimens", "primers", "runs"])
    input: "metadata.yml"
    run: R("""
            devtools::load_all("igseq")
            update_metadata_via_yaml("{input}", ".")
            """)
