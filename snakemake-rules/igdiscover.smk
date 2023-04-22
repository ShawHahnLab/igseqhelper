"""
Rules for running IgDiscover, starting from demultiplexed, adapter-trimmed,
merged IgM+ reads on a per-subject per-amplicon-type basis.

These rules are configured to allow different reference databases.

The default targets here will use SONAR's local copy of the Ramesh et al.
database (https://doi.org/10.3389/fimmu.2017.01407) from within our igseq
package.
"""

rule igdiscover_db_sonarramesh:
    # As per the IgDiscover manual: "The directory must contain the three files
    # V.fasta, D.fasta, J.fasta. These files contain the V, D, J gene sequences,
    # respectively. Even if you have only light chains in your data, a D.fasta
    # file needs to be provided; just use one with the heavy chain D gene
    # sequences."
    # So, we'll always include the IGHD sequences, plus whatever locus we're
    # working with here.  (And if that's IGH, it'll still just write the D
    # sequences once.)
    output: expand("analysis/igdiscover/sonarramesh/{{chain_type}}/{segment}.fasta", segment=["V", "D", "J"])
    params:
        outdir="analysis/igdiscover/sonarramesh/{chain_type}",
        locus=lambda w: {"mu": "IGH", "kappa": "IGK", "lambda": "IGL"}[w.chain_type]
    shell:
        """
            igseq vdj-gather sonarramesh/IGH/IGHD sonarramesh/{params.locus} -o {params.outdir}
        """

# KIMDB is only heavy chain so it's simpler to handle here than sonarramesh
# above
rule igdiscover_db_kimdb:
    output: expand("analysis/igdiscover/kimdb/mu/{segment}.fasta", segment=["V", "D", "J"])
    params: outdir="analysis/igdiscover/kimdb/mu"
    shell: "igseq vdj-gather kimdb -o {params.outdir}"

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
    conda: str(BASEDIR/"conda/igdiscover.yml")
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
    conda: str(BASEDIR/"conda/igdiscover.yml")
    threads: 20
    shell:
        """
            cd $(dirname {output.stats})/.. && igdiscover run --cores {threads}
        """
