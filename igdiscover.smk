"""
Rules for running IgDiscover, starting from demultiplexed, adapter-trimmed,
merged IgM+ reads on a per-subject per-amplicon-type basis.

These rules are configured to allow different reference databases.

The default targets here will use KIMDB originally published in Bernat et al.
2021 for heavy chain and SONAR's local copy of the Ramesh et al.  database
(https://doi.org/10.3389/fimmu.2017.01407) for light chain, both from within
our igseq package.
"""

import csv
import json
from datetime import (datetime, timezone)

# We'll use the standard outputs for V (and D) but will take the more
# stringently-filtered J output
# D will be ignored for light chain but is always there (which makes the
# snakemake rules here easy)
def input_for_igdiscover_as_germline(w):
    chain_type = {"IGH": "mu", "IGK": "kappa", "IGL": "lambda"}[w.locus]
    ref = "kimdb" if w.locus == "IGH" else "sonarramesh"
    targets = {
        "V": f"analysis/igdiscover/{ref}/{chain_type}/{w.subject}/final/database/V.fasta",
        "D": f"analysis/igdiscover/{ref}/{chain_type}/{w.subject}/final/database/D.fasta",
        "J": f"analysis/igdiscover/{ref}/{chain_type}/{w.subject}/custom_j_discovery/J.fasta",
        }
    return targets

rule igdiscover_as_germline:
    """Copy over IgDiscover output as the canonical germline reference for a subject"""
    output:
        fastas=expand("analysis/germline/{{subject}}.{{locus}}/{segment}.fasta", segment=["V", "D", "J"]),
        info="analysis/germline/{subject}.{locus}/info.csv"
    input: unpack(input_for_igdiscover_as_germline)
    run:
        shell("cp {input.V} {output.fastas[0]}")
        shell("cp {input.D} {output.fastas[1]}")
        shell("cp {input.J} {output.fastas[2]}")
        # Also, some basic information to help keep track of things
        dt_now = datetime.now(tz=timezone.utc).isoformat()
        # timestamp from IgDiscover log txt
        dt_log = "???"
        igdisc_log = Path(input.V).parent/"../../log.txt"
        if igdisc_log.exists():
            dt_log = datetime.fromtimestamp(
                Path(igdisc_log).stat().st_ctime, tz=timezone.utc).isoformat()
        # verison from stats.json
        igdisc_version = "???"
        igdisc_json = Path(input.V).parent/"../../stats/stats.json"
        if igdisc_json.exists():
            with open(igdisc_json) as f_in:
                igdisc_version = json.load(f_in)["version"]
        with open(output.info, "w") as f_out:
            writer = csv.writer(f_out, lineterminator="\n")
            writer.writerow(["key", "value"])
            writer.writerow(["time_run", dt_log])
            writer.writerow(["time_copied", dt_now])
            writer.writerow(["igdiscover_version", igdisc_version])
            writer.writerow(["V", input.V])
            writer.writerow(["D", input.D])
            writer.writerow(["J", input.J])


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
            Path(output[0]).symlink_to(Path(Path(input[0]).name)/inputs[0].name)
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
    conda: "igdiscover.yml"
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
    log:
        conda="analysis/igdiscover/{ref}/{chain_type}/{subject}/igdiscover_run.conda_build.txt"
    conda: "igdiscover.yml"
    threads: 20
    shell:
        """
            conda list --explicit > {log.conda}
            cd $(dirname {output.stats})/.. && igdiscover run --cores {threads}
        """

# An extra step for kappa J specifically:
#
# We sometimes see an implausible number of slight variations on IGKJ2 reported
# as novel alleles.  One IgDiscover author Martin Corcoran suggested more
# stringent options for the J discovery setting, specifically 100% for J
# coverage and 0.3 for the allele ratio (and suggested this may be related to
# the relatively low expression of specific genes such as IGKJ2).  In my
# testing it's the coverage setting that seems to do the trick here.
#
# from igdiscover discoverjd --help:
#
#  --j-coverage PERCENT  Require that the sequence covers at least PERCENT of
#                        the J gene. Default: 90 when --gene=J; 0 otherwise
#  --allele-ratio RATIO  Required allele ratio. Works only for genes named
#                        "NAME*ALLELE". Default: 0.2
rule custom_j_discovery:
    output:
        tab="analysis/igdiscover/{ref}/{chain_type}/{subject}/custom_j_discovery/J.tab",
        fasta="analysis/igdiscover/{ref}/{chain_type}/{subject}/custom_j_discovery/J.fasta"
    input:
        # (This command actually uses iteration-01's J.fasta and
        # filtered.tsv.gz but I don't have those listed as part of the
        # input/output paths in all this.  So I'll just request the output of
        # my igdiscover rule.)
        db_j="analysis/igdiscover/{ref}/{chain_type}/{subject}/final/database/J.fasta",
    params:
        db_j="analysis/igdiscover/{ref}/{chain_type}/{subject}/iteration-01/database/J.fasta",
        tab="analysis/igdiscover/{ref}/{chain_type}/{subject}/iteration-01/filtered.tsv.gz",
        jcov=100,
        ratio=0.3
    log:
        conda="analysis/igdiscover/{ref}/{chain_type}/{subject}/custom_j_discovery/conda_build.txt"
    conda: "igdiscover.yml"
    shell:
        """
            arg_j_cov=""
            if [[ "{params.jcov}" != "" ]]; then
                arg_j_cov="--j-coverage {params.jcov}"
            fi
            arg_allele_ratio=""
            if [[ "{params.ratio}" != "" ]]; then
                arg_allele_ratio="--allele-ratio {params.ratio}"
            fi
            conda list --explicit > {log.conda}
            igdiscover discoverjd --database {params.db_j} --gene J $arg_j_cov $arg_allele_ratio {params.tab} --fasta {output.fasta} > {output.tab}
        """

rule upstream:
    """Derive a consensus sequence for the upstream portion (UTR+Leader) of each V gene."""
    output: "analysis/igdiscover-upstream/{ref}.{chain_type}.{subject}/upstream.fasta"
    input: "analysis/igdiscover/{ref}/{chain_type}/{subject}/final/database/V.fasta"
    params:
        table="analysis/igdiscover/{ref}/{chain_type}/{subject}/final/filtered.tsv.gz"
    conda: "gkhlab-igdiscover22.yml"
    log:
        main="analysis/igdiscover-upstream/{ref}.{chain_type}.{subject}/log.txt",
        conda="analysis/igdiscover-upstream/{ref}.{chain_type}.{subject}/conda_build.txt"
    shell:
        """
            conda list --explicit > {log.conda}
            igdiscover upstream {params.table} > {output} 2> {log.main}
        """
