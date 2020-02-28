"""
Rules for running IgDisover, starting from demultiplexed, adapter-trimmed reads
on a per-specimen per-amplicon-type basis.  IgDiscover does its own
read-merging using PEAR so we won't use pRESTO's merge reads here.

The databases will be prepared from SONAR's local copy of the Ramesh et al.
database (https://doi.org/10.3389/fimmu.2017.01407).
"""

from igseq.igdiscover import CHAINS, CHAIN_TYPES, SEGMENTS

# We'll run IgDiscover on all IgM+ specimens, so, some combination of heavy.mu,
# light.lambda, light.kappa.
TARGET_IGDISCOVER_INIT = amplicon_files(
    "igdiscover/{chain}.{chain_type}/{specimen}/igdiscover.yaml", SAMPLES, "IgM+")
TARGET_IGDISCOVER_ALL = amplicon_files(
    "igdiscover/{chain}.{chain_type}/{specimen}/stats/stats.json", SAMPLES, "IgM+")

rule all_igdiscover_init:
    input: TARGET_IGDISCOVER_INIT

rule all_igdiscover:
    input: TARGET_IGDISCOVER_ALL

def input_igdiscover_db(w):
    """Gather VDJ germline DB from SONAR's BU_DD files.

    As per the IgDiscover manual: "The directory must contain the three files
    V.fasta, D.fasta, J.fasta. These files contain the V, D, J gene sequences,
    respectively. Even if you have only light chains in your data, a D.fasta
    file needs to be provided; just use one with the heavy chain D gene
    sequences."
    """
    # special case: IgDiscover wants a dummy D for a light chain so we'll use
    # H's D.
    chain_type = CHAIN_TYPES[w.chain_type]
    if CHAIN_TYPES[w.chain_type] != "H" and w.segment == "D":
        chain_type = "H"
    return expand("SONAR/germDB/Ig{x}{z}_BU_DD.fasta", x=chain_type, z=w.segment)

rule igdiscover_db:
    output: "igdiscover/{chain}.{chain_type}/{segment}.fasta"
    input: input_igdiscover_db
    shell:
        """
            cp {input} {output}
        """

# SNAKEMAKES ON SNAKEMAKES
rule igdiscover_init:
    output: "igdiscover/{chain}.{chain_type}/{specimen}/igdiscover.yaml"
    input:
        # IgDiscover always wants a "D" even for light
        db_v="igdiscover/{chain}.{chain_type}/V.fasta",
        db_d="igdiscover/{chain}.{chain_type}/D.fasta",
        db_j="igdiscover/{chain}.{chain_type}/J.fasta",
        r1="presto/data/{chain}.{chain_type}/{specimen}.R1.fastq",
        r2="presto/data/{chain}.{chain_type}/{specimen}.R2.fastq"
    params: iterations=5
    shell:
        """
            rmdir $(dirname {output})
            igdiscover init \
                --db $(dirname {input.db_v}) \
                --reads {input.r1} \
                $(dirname {output})
            sed -i 's/^iterations: 1$/iterations: {params.iterations}/' {output}
        """

rule igdiscover_run:
    output: "igdiscover/{chain}.{chain_type}/{specimen}/stats/stats.json"
    input: "igdiscover/{chain}.{chain_type}/{specimen}/igdiscover.yaml"
    threads: 20
    shell:
        """
            cd $(dirname {output})/.. && igdiscover run --cores {threads}
        """
