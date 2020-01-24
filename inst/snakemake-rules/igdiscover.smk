"""
Rules for running IgDisover, starting from demultiplexed, adapter-trimmed reads.

The databases will be prepared from the IMGT reference for Rhesus Macaque in
the igseq directory.
"""

CHAINS = {"heavy": "H", "light": "L"}
SEGMENTS = {"heavy": ["V", "D", "J"], "light": ["V", "J"]}

rule igdiscover_db:
    output: "igdiscover/{chain}/{segment}.fasta"
    input: lambda w: expand("igseq/inst/reference/imgt/IG{chain}{segment}.fasta", chain = CHAINS[w.chain], segment = w.segment)
    shell:
        """
            awk '{{if(NR%2==0){{gsub("\\\\.","");print $0}}else{{print $0}}}}' < {input} > {output}
        """

# SNAKEMAKES ON SNAKEMAKES
rule igdiscover_init:
    output: "igdiscover/{chain}/{sample}/igdiscover.yaml"
    input:
        db=lambda w: expand("igdiscover/{chain}/{segment}.fasta", chain = w.chain, segment = SEGMENTS[w.chain]),
        r1="presto/data/{sample}_1.fastq"
    shell:
        """
            rmdir $(dirname {output})
            igdiscover init \
                --db $(dirname {input.db[0]}) \
                --reads {input.r1} \
                $(dirname {output})
        """

rule igdiscover_run:
    output: "igdiscover/{chain}/{sample}/stats/stats.json"
    input: "igdiscover/{chain}/{sample}/igdiscover.yaml"
    threads: 20
    shell:
        """
            cd $(dirname {output})/.. && igdiscover run --cores {threads}
        """
