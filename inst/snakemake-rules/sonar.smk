import igseq.sonar

TARGET_SONAR_PREP = \
    expand("sonar-analysis/heavy/{sample}/{sample}.fastq", sample = SAMPLES_H) + \
    expand("sonar-analysis/light/{sample}/{sample}.fastq", sample = SAMPLES_L)

TARGET_SONAR_MODULE_1 = \
    expand("sonar-analysis/heavy/{sample}/output/sequences/nucleotide/{sample}_goodVJ_unique.fa", sample = SAMPLES_H) + \
    expand("sonar-analysis/light/{sample}/output/sequences/nucleotide/{sample}_goodVJ_unique.fa", sample = SAMPLES_L)

TARGET_SONAR_MODULE_2_AUTO = \
    expand("sonar-analysis/heavy/{sample}/output/tables/{sample}_lineages.txt", sample = SAMPLES_H) + \
    expand("sonar-analysis/light/{sample}/output/tables/{sample}_lineages.txt", sample = SAMPLES_L)

TARGET_SONAR_MODULE_2_ID_DIV = \
    expand("sonar-analysis/heavy/{sample}/output/tables/{sample}_goodVJ_unique_id-div.tab", sample = SAMPLES_H) + \
    expand("sonar-analysis/light/{sample}/output/tables/{sample}_goodVJ_unique_id-div.tab", sample = SAMPLES_L)

TARGET_SONAR_MODULE_2_MANUAL = \
    expand("sonar-analysis/heavy/{sample}/output/sequences/nucleotide/{sample}_islandSeqs.fa", sample = SAMPLES_H) + \
    expand("sonar-analysis/light/{sample}/output/sequences/nucleotide/{sample}_islandSeqs.fa", sample = SAMPLES_L)

TARGET_SONAR_MODULE_3 = \
    expand("sonar-analysis/{chain}/longitudinal/output/longitudinal_igphyml.tree", chain = ["heavy", "light"])

rule all_sonar_prep_input:
    input: TARGET_SONAR_PREP

rule all_sonar_module_1:
    input: TARGET_SONAR_MODULE_1

rule all_sonar_module_2_auto:
    input: TARGET_SONAR_MODULE_2_AUTO

rule all_sonar_module_2_id_div:
    input: TARGET_SONAR_MODULE_2_ID_DIV

rule all_sonar_module_2_manual:
    input: TARGET_SONAR_MODULE_2_MANUAL

rule all_sonar_module_3:
    input: TARGET_SONAR_MODULE_3

ruleorder: sonar_prep_input_from_presto > sonar_prep_input

# Concatenate all reads found across runs for a single sample into one
# fastq, inside a subdirectory for that sample's SONAR analysis.
rule sonar_prep_input:
    output: "sonar-analysis/{chain}/{sample}/{sample}.fastq"
    input: expand("5-primertrim/{run}/{{sample}}.fastq.gz", run = RUNS.keys())
    shell: "zcat {input} > {output}"

# New version!  use pRESTO instead of our custom prep steps
rule sonar_prep_input_from_presto:
    output: "sonar-analysis/{chain}/{sample}/{sample}.fastq"
    input: "presto/qual/{sample}-FWD_primers-pass.fastq"
    shell: "ln {input} {output}"

# Prepare germline references for SONAR.  Note that it expects a particular
# sequence naming scheme for genes and module 2 will fail if we don't follow
# it.
rule sonar_gather_germline:
    output: "sonar-analysis/{chain}/germline.{segment}.fasta"
    input:
        imgt=lambda w: str(RINST/"reference/imgt/IG%s%s.fasta") % ({"heavy": "H", "light": "L"}[w.chain], w.segment),
        extra_v=str(RINST/"reference/10.1016_j.cell.2019.06.030/tables3c.{chain}.fasta")
    run:
        igseq.sonar.gather_germline(input.imgt, input.extra_v, output[0], wildcards.segment)

rule sonar_gather_mature:
    output: "sonar-analysis/{chain}/mab.fasta"
    input: str(RINST/"reference/{chain}.mature.fasta")
    shell: "cp {input} {output}"

# use --jmotif arguent to sonar finalize?
# ("Conserved nucleotide sequence indicating the start of FWR4 on the J gene.
# Defaults to either TGGGG for heavy chains or TT[C|T][G|A]G for light chains;
# set manually for custom light chain libraries or for species with a different
# motif.")
# In the files from IMGT, we do have TGGGG in all of the IGHJ.fasta sequences
# and TTCGG in all of the IGLJ.fasta sequences, so maybe that's OK as the
# default? Also, our amplicons seem to end with just a tiny fragment of this
# conserved region anyway.
rule sonar_module_1:
    output: "sonar-analysis/{chain}/{sample}/output/sequences/nucleotide/{sample}_goodVJ_unique.fa"
    input:
        fq="sonar-analysis/{chain}/{sample}/{sample}.fastq",
        ref_hv="sonar-analysis/heavy/germline.V.fasta",
        ref_hd="sonar-analysis/heavy/germline.D.fasta",
        ref_hj="sonar-analysis/heavy/germline.J.fasta",
        ref_lv="sonar-analysis/light/germline.V.fasta",
        ref_lj="sonar-analysis/light/germline.J.fasta"
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        cluster_id_fract=.97,
        cluster_min2=2
    shell:
        """
            cd $(dirname {input.fq})
            if [[ {wildcards.chain} == "heavy" ]]; then
                true
                libv_arg="--lib ../$(basename {input.ref_hv})"
                libd_arg="--dlib ../$(basename {input.ref_hd})"
                libj_arg="--jlib ../$(basename {input.ref_hj})"
            else
                libv_arg="--lib ../$(basename {input.ref_lv})"
                libd_arg="--noD"
                libj_arg="--jlib ../$(basename {input.ref_lj})"
            fi
            sonar blast_V --fasta ../../../{input.fq} $libv_arg  --derep --threads {threads}
            sonar blast_J $libd_arg $libj_arg --noC --threads {threads}
            sonar finalize --noclean --threads {threads}
            sonar cluster_sequences --id {params.cluster_id_fract} --min2 {params.cluster_min2}
        """

### Non-interactive aspect of module 2

# Automated intradonor analysis and grouping
rule sonar_module_2:
    output: "sonar-analysis/{chain}/{sample}/output/tables/{sample}_lineages.txt"
    input:
        module1="sonar-analysis/{chain}/{sample}/output/sequences/nucleotide/{sample}_goodVJ_unique.fa",
        mab="sonar-analysis/{chain}/mab.fasta",
        ref_hv="sonar-analysis/heavy/germline.V.fasta",
        ref_hd="sonar-analysis/heavy/germline.D.fasta",
        ref_hj="sonar-analysis/heavy/germline.J.fasta",
        ref_lv="sonar-analysis/light/germline.V.fasta",
        ref_lj="sonar-analysis/light/germline.J.fasta",
        ref_mab="sonar-analysis/{chain}/mab.fasta",
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        vh_id=ALLELES["heavy"]["V"],
        vl_id=ALLELES["light"]["V"],
        jh_id=ALLELES["heavy"]["J"],
        jl_id=ALLELES["light"]["J"]
    shell:
        """
            cd $(dirname {input.module1})/../../..
            libv=../germline.V.fasta
            mab=../mab.fasta
            if [[ {wildcards.chain} == heavy ]]; then
                vid="{params.vh_id}"
                jid="{params.jh_id}"
            else
                vid="{params.vl_id}"
                jid="{params.jl_id}"
            fi
            # 2) Automated intradonor analysis
            # Note though what the paper says about the intradonor results
            # overlaid on the id-div results in figure 3: "Two thirds of these
            # transcripts are found in the high-identity island; the remaining
            # third in the main body of transcripts at ~70% identity are false
            # positives. This is a typical result, showing why multiple tools
            # for lineage determination are included in SONAR and manual
            # curation is strongly advised."
            sonar intradonor --n "$mab" --v "$vid" --lib "$libv" --threads {threads}
            # 3) Grouping by V and J and clustering by CDR3
            # We could use --natives as well but I'm not 100% sure what it
            # expects in the FASTA in that case.  Just CDR3?  If so we could
            # chop out that segment from Ryan's mature references and provide
            # it here.
            sonar groups -v "$vid" -j "$jid" -t {threads}
        """

### Manually-guided aspect of Module 2
# Identity/Divergence plots and island identification
# The first part-- calculating the identity/divergence coordinates-- is
# automatic.  The second part-- selecting points on the scatterplots-- is
# manual and requires X11. (The third is just extracting the selected sequences
# from the original FASTA.)

# Part 1 of 3: Calculate % identity to each of the mature antibodies.
rule sonar_module_2_id_div:
    output:
        iddiv="sonar-analysis/{chain}/{sample}/output/tables/{sample}_goodVJ_unique_id-div.tab"
    input:
        fasta="sonar-analysis/{chain}/{sample}/output/sequences/nucleotide/{sample}_goodVJ_unique.fa",
        mab="sonar-analysis/{chain}/mab.fasta",
        germline_v="sonar-analysis/{chain}/germline.V.fasta"
    singularity: "docker://scharch/sonar"
    threads: 4
    shell:
        """
            cd $(dirname {input.fasta})/../../..
            libv=../../../{input.germline_v}
            mab=../../../{input.mab}
            sonar id-div -g "$libv" -a "$mab" -t {threads}
        """

# Part 2 of 3: Select specific clusters from id/div plots
# NOTE this step is interactive over X11
# haven't been able to get this part working with snakemake's
# singularity/docker stuff
rule sonar_module_2_id_div_island:
    output:
        seqids="sonar-analysis/{chain}/{sample}/output/tables/islandSeqs.txt"
    input:
        iddiv="sonar-analysis/{chain}/{sample}/output/tables/{sample}_goodVJ_unique_id-div.tab"
    params:
        mab="=a", # =a means use all antibodies in mab input file
    shell:
        """
            set +e
            set +u
            source ~/miniconda3/bin/activate sonar
            cd $(dirname {input.tab})/../..
            sonar get_island ../../../{input.iddiv} --mab "{params.mab}"
        """

# Part 3 of 3: Extract those goodVJ cluster sequences given in the id/div lists
rule sonar_module_2_id_div_getfasta:
    output:
        fasta="sonar-analysis/{chain}/{sample}/output/sequences/nucleotide/{sample}_islandSeqs.fa"
    input:
        fasta="sonar-analysis/{chain}/{sample}/output/sequences/nucleotide/{sample}_goodVJ_unique.fa",
        seqids="sonar-analysis/{chain}/{sample}/output/tables/islandSeqs.txt"
    singularity: "docker://scharch/sonar"
    shell:
        """
            cd $(dirname {input.fasta})/../../..
            # Pull out sequences using our selected island from the above
            # id-div step.  I think this is essentially "seqmagick convert
            # --include-from-file ..."
            sonar getfasta \
                -l ../../../{input.seqids} \
                -f ../../../{input.fasta} \
                -o ../../../{output.fasta}
        """


rule sonar_module_3_collect_light:
    output:
        collected="sonar-analysis/light/longitudinal/output/sequences/nucleotide/longitudinal-collected.fa"
    input:
        WK12L="sonar-analysis/light/WK12L/output/sequences/nucleotide/WK12L_islandSeqs.fa",
        WK16L="sonar-analysis/light/WK16L/output/sequences/nucleotide/WK16L_islandSeqs.fa",
        WK24L="sonar-analysis/light/WK24L/output/sequences/nucleotide/WK24L_islandSeqs.fa",
        WK48L="sonar-analysis/light/WK48L/output/sequences/nucleotide/WK48L_islandSeqs.fa",
        WK56L="sonar-analysis/light/WK56L/output/sequences/nucleotide/WK56L_islandSeqs.fa",
        WK65L="sonar-analysis/light/WK65L/output/sequences/nucleotide/WK65L_islandSeqs.fa"
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd="sonar-analysis/light/longitudinal"
    shell:
        """
            cd {params.wd}
            sonar merge_time \
                --seqs ../../../{input.WK12L} --labels WK12L \
                --seqs ../../../{input.WK16L} --labels WK16L \
                --seqs ../../../{input.WK24L} --labels WK24L \
                --seqs ../../../{input.WK48L} --labels WK48L \
                --seqs ../../../{input.WK56L} --labels WK56L \
                --seqs ../../../{input.WK65L} --labels WK65L
        """

rule sonar_module_3_igphyml_light:
    output:
        tree="sonar-analysis/light/longitudinal/output/longitudinal_igphyml.tree",
        inferred_nucl="sonar-analysis/light/longitudinal/output/sequences/nucleotide/longitudinal_inferredAncestors.fa",
        inferred_prot="sonar-analysis/light/longitudinal/output/sequences/amino_acid/longitudinal_inferredAncestors.fa",
        stats="sonar-analysis/light/longitudinal/output/logs/longitudinal_igphyml_stats.txt"
    input:
        collected="sonar-analysis/light/longitudinal/output/sequences/nucleotide/longitudinal-collected.fa",
        ref_hv="sonar-analysis/light/germline.V.fasta",
        natives=str(RINST/"reference/light.mature.fasta")
    #singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd="sonar-analysis/light/longitudinal",
        vl_id=ALLELES["light"]["V"],
        args="-f"
    shell:
        """
            set +e
            set +u
            source ~/miniconda3/bin/activate sonar
            cd {params.wd}
            sonar igphyml \
                -v '{params.vl_id}' \
                --lib ../germline.V.fasta \
                --natives {input.natives} \
                {params.args}
        """

rule sonar_module_3_collect_heavy:
    output:
        collected="sonar-analysis/heavy/longitudinal/output/sequences/nucleotide/longitudinal-collected.fa"
    input:
        WK12H="sonar-analysis/heavy/WK12H/output/sequences/nucleotide/WK12H_islandSeqs.fa",
        WK16H="sonar-analysis/heavy/WK16H/output/sequences/nucleotide/WK16H_islandSeqs.fa",
        WK16H2="sonar-analysis/heavy/WK16H2/output/sequences/nucleotide/WK16H2_islandSeqs.fa",
        WK24H="sonar-analysis/heavy/WK24H/output/sequences/nucleotide/WK24H_islandSeqs.fa",
        WK48H="sonar-analysis/heavy/WK48H/output/sequences/nucleotide/WK48H_islandSeqs.fa",
        WK56H="sonar-analysis/heavy/WK56H/output/sequences/nucleotide/WK56H_islandSeqs.fa",
        WK56H2="sonar-analysis/heavy/WK56H2/output/sequences/nucleotide/WK56H2_islandSeqs.fa",
        WK65H="sonar-analysis/heavy/WK65H/output/sequences/nucleotide/WK65H_islandSeqs.fa"
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd="sonar-analysis/heavy/longitudinal"
    shell:
        """
            # Note that here we have a couple of replicates (WK16H2, WK56H2)
            # that we're putting under the same timepoint label as the first of
            # each pair.
            cd {params.wd}
            sonar merge_time \
                --seqs ../../../{input.WK12H} --labels WK12H \
                --seqs ../../../{input.WK16H} --labels WK16H \
                --seqs ../../../{input.WK16H2} --labels WK16H \
                --seqs ../../../{input.WK24H} --labels WK24H \
                --seqs ../../../{input.WK48H} --labels WK48H \
                --seqs ../../../{input.WK56H} --labels WK56H \
                --seqs ../../../{input.WK56H2} --labels WK56H \
                --seqs ../../../{input.WK65H} --labels WK65H
        """

rule sonar_module_3_igphyml_heavy:
    output:
        tree="sonar-analysis/heavy/longitudinal/output/longitudinal_igphyml.tree",
        inferred_nucl="sonar-analysis/heavy/longitudinal/output/sequences/nucleotide/longitudinal_inferredAncestors.fa",
        inferred_prot="sonar-analysis/heavy/longitudinal/output/sequences/amino_acid/longitudinal_inferredAncestors.fa",
        stats="sonar-analysis/heavy/longitudinal/output/logs/longitudinal_igphyml_stats.txt"
    input:
        collected="sonar-analysis/heavy/longitudinal/output/sequences/nucleotide/longitudinal-collected.fa",
        ref_hv="sonar-analysis/heavy/germline.V.fasta",
        natives=str(RINST/"reference/heavy.mature.fasta")
    #singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd="sonar-analysis/heavy/longitudinal",
        vh_id=ALLELES["heavy"]["V"],
        args="-f"
    shell:
        """
            set +e
            set +u
            source ~/miniconda3/bin/activate sonar
            cd {params.wd}
            sonar igphyml \
                -v '{params.vh_id}' \
                --lib ../germline.V.fasta \
                --natives {input.natives} \
                {params.args}
        """
