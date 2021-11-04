"""
Rules for running IgDisover, starting from demultiplexed, adapter-trimmed,
merged IgM+ reads on a per-subject per-amplicon-type basis.

These rules are configured to allow different reference databases.

The default targets here will use SONAR's local copy of the Ramesh et al.
database (https://doi.org/10.3389/fimmu.2017.01407).  This requires the SONAR
repository in the working directory.
"""

def input_igdiscover_db_sonar(w):
    """Gather VDJ germline DB from SONAR's BU_DD files.

    As per the IgDiscover manual: "The directory must contain the three files
    V.fasta, D.fasta, J.fasta. These files contain the V, D, J gene sequences,
    respectively. Even if you have only light chains in your data, a D.fasta
    file needs to be provided; just use one with the heavy chain D gene
    sequences."
    """
    # Mapping of chain types to heavy/light short specifiers per locus.
    chain_types = {
        "alpha": "H",
        "delta": "H",
        "gamma": "H",
        "mu": "H",
        "epsilon": "H",
        "kappa": "K",
        "lambda": "L"}
    # special case: IgDiscover wants a dummy D for a light chain so we'll use
    # H's D.
    chain_type = chain_types[w.chain_type]
    if chain_types[w.chain_type] != "H" and w.segment == "D":
        chain_type = "H"
    return expand("SONAR/germDB/Ig{x}{z}_BU_DD.fasta", x=chain_type, z=w.segment)

rule igdiscover_db_sonarramesh:
    output: "analysis/igdiscover/SONARRamesh/{chain_type}/{segment}.fasta"
    input: input_igdiscover_db_sonar
    shell:
        """
            cp {input} {output}
        """

rule catsubjects:
    output: "analysis/samples-by-subject/igm/{subject}.{chain_type}.fastq.gz"
    input: "analysis/samples-by-subject/igm/{subject}.{chain_type}"
    run:
        # usually there will be just one file destined for IgDiscover for a
        # given subject, but once in an a while we have more than one.  We'll
        # symlink when there's only one and concatenate otherwise.
        inputs = list(Path(input[0]).glob("*.fastq.gz"))
        if len(inputs) == 1:
            Path(output[0]).resolve().symlink_to(Path(inputs[0]).resolve())
        elif inputs:
            shell("cat {input}/*.fastq.gz > {output}")
        else:
            raise ValueError

# SNAKEMAKES ON SNAKEMAKES
rule igdiscover_init:
    output:
        yaml="analysis/igdiscover/{ref}/{chain_type}/{subject}/igdiscover.yaml",
        reads="analysis/igdiscover/{ref}/{chain_type}/{subject}/reads.fastq.gz",
    input:
        # IgDiscover always wants a "D" even for light
        db_v="analysis/igdiscover/{ref}/{chain_type}/V.fasta",
        db_d="analysis/igdiscover/{ref}/{chain_type}/D.fasta",
        db_j="analysis/igdiscover/{ref}/{chain_type}/J.fasta",
        reads="analysis/samples-by-subject/igm/{subject}.{chain_type}.fastq.gz",
    conda: str(BASEDIR/"igdiscover.yml")
    params:
        stranded="true",
        iterations=5
    shell:
        """
            rmdir $(dirname {output.yaml})
            igdiscover init \
                --db $(dirname {input.db_v}) \
                --single-reads {input.reads} \
                $(dirname {output.yaml})
            # This is a sloppy way of handling the YAML!  Re-work this probably
            # by parsing in Python and writing back out.
            sed -i 's/^iterations: 1$/iterations: {params.iterations}/' {output.yaml}
            sed -i 's/^stranded: false$/stranded: {params.stranded}/' {output.yaml}
        """

rule igdiscover_run:
    output:
        stats="analysis/igdiscover/{ref}/{chain_type}/{subject}/stats/stats.json",
        db_v="analysis/igdiscover/{ref}/{chain_type}/{subject}/final/database/V.fasta",
        db_d="analysis/igdiscover/{ref}/{chain_type}/{subject}/final/database/D.fasta",
        db_j="analysis/igdiscover/{ref}/{chain_type}/{subject}/final/database/J.fasta"
    input:
        yaml="analysis/igdiscover/{ref}/{chain_type}/{subject}/igdiscover.yaml",
        r1="analysis/igdiscover/{ref}/{chain_type}/{subject}/reads.fastq.gz",
    conda: str(BASEDIR/"igdiscover.yml")
    threads: 20
    shell:
        """
            cd $(dirname {output.stats})/.. && igdiscover run --cores {threads}
        """
