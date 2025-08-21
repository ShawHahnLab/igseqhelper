"""
Rules for reporting on IgDiscover output.
"""

# See also the final/dendrogram_{segment}.pdf files from IgDiscover
rule report_igdiscover_tree:
    output: "analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.nex"
    input:
        after="analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.fasta",
        before="analysis/igdiscover/{ref}/{chain_type}/{segment}.fasta"
    log:
        conda="analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.nex.conda_build.txt"
    conda: "envs/igseq.yaml"
    conda: "envs/igseq.yaml"
    shell:
        """
            conda list --explicit > {log.conda}
            igseq tree before={input.before} after={input.after} {output}
        """

rule report_igdiscover_final_db:
    output: "analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.fasta"
    input: "analysis/igdiscover/{ref}/{chain_type}/{subject}/final/database/{segment}.fasta"
    shell: "cp {input} {output}"
