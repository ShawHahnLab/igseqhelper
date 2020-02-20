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
    # special case: IgDiscover wants an empty D for the light chain
    if CHAIN_TYPES[w.chain_type] != "H" and w.segment == "D":
        return "/dev/null"
    return expand("SONAR/germDB/Ig{x}{z}_BU_DD.fasta", x=CHAIN_TYPES[w.chain_type], z=w.segment)

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
        db=expand("igdiscover/{{chain}}.{{chain_type}}/{segment}.fasta", segment=SEGMENTS["heavy"]),
        r1="presto/data/{chain}.{chain_type}/{specimen}.R1.fastq",
        r2="presto/data/{chain}.{chain_type}/{specimen}.R2.fastq"
    shell:
        """
            rmdir $(dirname {output})
            igdiscover init \
                --db $(dirname {input.db[0]}) \
                --reads {input.r1} \
                $(dirname {output})
        """

rule igdiscover_run:
    output: "igdiscover/{chain}.{chain_type}/{specimen}/stats/stats.json"
    input: "igdiscover/{chain}.{chain_type}/{specimen}/igdiscover.yaml"
    threads: 20
    shell:
        """
            cd $(dirname {output})/.. && igdiscover run --cores {threads}
        """
