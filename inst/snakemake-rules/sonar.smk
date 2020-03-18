import igseq.sonar
from igseq.data import MetadataError

IGG_CHAINS = [entry["Chain"] for entry in SAMPLES.values() if "IgG+" in entry["SpecimenAttrs"]["CellType"]]
IGG_CHAINTYPES = [entry["Type"] for entry in SAMPLES.values() if "IgG+" in entry["SpecimenAttrs"]["CellType"]]
IGG_SUBJECTS = [entry["SpecimenAttrs"]["Subject"] for entry in SAMPLES.values() if "IgG+" in entry["SpecimenAttrs"]["CellType"]]
IGG_SPECIMENS = [entry["Specimen"] for entry in SAMPLES.values() if "IgG+" in entry["SpecimenAttrs"]["CellType"]]

TARGET_SONAR_PREP = expand(
    "analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/{specimen}.fastq",
    zip, subject=IGG_SUBJECTS, chain=IGG_CHAINS, chain_type=IGG_CHAINTYPES, specimen=IGG_SPECIMENS)
TARGET_SONAR_MODULE_1 = expand(
    "analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
    zip, subject=IGG_SUBJECTS, chain=IGG_CHAINS, chain_type=IGG_CHAINTYPES, specimen=IGG_SPECIMENS)
TARGET_SONAR_MODULE_2_AUTO = expand(
    "analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/tables/{specimen}_lineages.txt",
    zip, subject=IGG_SUBJECTS, chain=IGG_CHAINS, chain_type=IGG_CHAINTYPES, specimen=IGG_SPECIMENS)
TARGET_SONAR_MODULE_2_ID_DIV = expand(
    "analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab",
    zip, subject=IGG_SUBJECTS, chain=IGG_CHAINS, chain_type=IGG_CHAINTYPES, specimen=IGG_SPECIMENS)
TARGET_SONAR_MODULE_2_MANUAL = expand(
    "analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/sequences/nucleotide/{specimen}_islandSeqs.fa",
    zip, subject=IGG_SUBJECTS, chain=IGG_CHAINS, chain_type=IGG_CHAINTYPES, specimen=IGG_SPECIMENS)
TARGET_SONAR_MODULE_3 = expand(
    "analysis/sonar/{subject}/{chain}.{chain_type}/longitudinal/{subject}/output/longitudinal_igphyml.tree",
    zip, chain=IGG_CHAINS, chain_type=IGG_CHAINTYPES, subject=IGG_SUBJECTS) 

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

def sonar_gather_germline_inputs(wildcards):
    """Get IgDiscover per-specimen paths for a given per-subject output."""
    # IgG+ will have gamma for heavy, but the germline will have come from mu.
    return igseq.sonar.igdiscover_final_db(
        SAMPLES, wildcards.subject, wildcards.chain, wildcards.chain_type, wildcards.segment)

rule sonar_prep_input_from_presto:
    """Gather input sequences from pRESTO's initial processing."""
    output: "analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/{specimen}.fastq"
    input: "analysis/presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq"
    shell: "ln {input} {output}"

rule sonar_gather_germline:
    """Gather germline alleles from IgDiscover's analysis of IgM+ specimens.

    Prepare germline references for SONAR.  Note that it expects a particular
    sequence naming scheme for genes and module 2 will fail if we don't follow
    it.  This is written to handle the possibility of multiple specimens in the
    IgDiscover input and pool them, though as of this writing we just use one.

    We'll need to take germline references from corresponding specimens for the
    same subject, using mu for heavy chain and lambda/kappa for light, though
    the corresponding chain types for the samples here for SONAR will be gamma
    for heavy and lambda/kappa for light.
    """
    output: "analysis/sonar/{subject}/{chain}.{chain_type}/germline.{segment}.fasta"
    input: sonar_gather_germline_inputs
    run:
        igseq.sonar.gather_germline(input, output[0])

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
    output:
        fasta="analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        rearr="analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/tables/{specimen}_rearrangements.tsv",
    input: unpack(igseq.sonar.sonar_module_1_inputs)
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        cluster_id_fract=.97,
        cluster_min2=2,
        # What we give for the V(D)J command-line argments depends on which
        # chain we're doing.  We could handle this during the rule itself with
        # a run: directive but then it won't let us use singularity.
        libv_arg=lambda w, input: "--lib ../" + str(Path(input.V).name),
        libd_arg=lambda w, input: "D" in input and "--dlib ../" + str(Path(input.D).name) or "--noD",
        libj_arg=lambda w, input: "--jlib ../" + str(Path(input.J).name)
    shell:
        """
            cd $(dirname {input.fastq})
            sonar blast_V --fasta ../../../../{input.fastq} {params.libv_arg} --derep --threads {threads}
            sonar blast_J {params.libd_arg} {params.libj_arg} --noC --threads {threads}
            sonar finalize --noclean --threads {threads}
            sonar cluster_sequences --id {params.cluster_id_fract} --min2 {params.cluster_min2}
        """

rule sonar_gather_mature:
    """Get heavy or light chain mature antibody sequences for a subject."""
    output: "analysis/sonar/{subject}/{chain}.{chain_type}/mab.fasta"
    run:
        igseq.sonar.gather_mature(
            ANTIBODY_ISOLATES,
            wildcards.subject,
            wildcards.chain, 
            output[0])

### Non-interactive aspect of module 2

def sonar_allele_v(wildcards):
    """Get the sequence ID for the expected germline V allele"""
    seqids = get_antibody_alleles(ANTIBODY_LINEAGES, wildcards.subject, wildcards.chain, "V")
    if len(seqids) != 1:
        raise MetadataError("More than one antibody lineage for subject: %s" % wildcards.subject)
    return seqids[0]

def sonar_allele_j(wildcards):
    """Get the sequence ID for the expected germline J allele"""
    seqids = get_antibody_alleles(ANTIBODY_LINEAGES, wildcards.subject, wildcards.chain, "J")
    if len(seqids) != 1:
        raise MetadataError("More than one antibody lineage for subject: %s" % wildcards.subject)
    return seqids[0]

# Automated intradonor analysis and grouping
rule sonar_module_2:
    output: "analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/tables/{specimen}_lineages.txt"
    input: unpack(igseq.sonar.sonar_module_2_inputs)
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        v_id=sonar_allele_v,
        j_id=sonar_allele_j
    shell:
        """
            # Start from the top-level directory for this specimen's analysis
            cd $(dirname {input.module1})/../../..
            libv=../germline.V.fasta
            mab=../mab.fasta
            # 2) Automated intradonor analysis
            # Note though what the paper says about the intradonor results
            # overlaid on the id-div results in figure 3: "Two thirds of these
            # transcripts are found in the high-identity island; the remaining
            # third in the main body of transcripts at ~70% identity are false
            # positives. This is a typical result, showing why multiple tools
            # for lineage determination are included in SONAR and manual
            # curation is strongly advised."
            sonar intradonor --n "$mab" --v "{params.v_id}" --lib "$libv" --threads {threads}
            # 3) Grouping by V and J and clustering by CDR3
            # We could use --natives as well but I'm not 100% sure what it
            # expects in the FASTA in that case.  Just CDR3?  If so we could
            # chop out that segment from Ryan's mature references and provide
            # it here.
            sonar groups -v "{params.v_id}" -j "{params.j_id}" -t {threads}
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
        iddiv="analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab"
    input:
        fasta="analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        mab="analysis/sonar/{subject}/{chain}.{chain_type}/mab.fasta",
        germline_v="analysis/sonar/{subject}/{chain}.{chain_type}/germline.V.fasta"
    singularity: "docker://scharch/sonar"
    threads: 4
    shell:
        """
            cd $(dirname {input.fasta})/../../..
            libv=../germline.V.fasta
            mab=../mab.fasta
            sonar id-div -g "$libv" -a "$mab" -t {threads}
        """

# Part 2 of 3: Select specific clusters from id/div plots
# NOTE this step is interactive over X11
# haven't been able to get this part working with snakemake's
# singularity/docker stuff
rule sonar_module_2_id_div_island:
    output:
        seqids="analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/tables/islandSeqs.txt"
    input:
        iddiv="analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab"
    params:
        mab="=a", # =a means use all antibodies in mab input file
    shell:
        """
            set +e
            set +u
            source ~/miniconda3/bin/activate sonar
            cd $(dirname {input.iddiv})/../..
            sonar get_island ../../../../{input.iddiv} --mab "{params.mab}"
        """

# Part 3 of 3: Extract those goodVJ cluster sequences given in the id/div lists
rule sonar_module_2_id_div_getfasta:
    output:
        fasta="analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/sequences/nucleotide/{specimen}_islandSeqs.fa"
    input:
        fasta="analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        seqids="analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/tables/islandSeqs.txt"
    singularity: "docker://scharch/sonar"
    shell:
        """
            cd $(dirname {input.fasta})/../../..
            # Pull out sequences using our selected island from the above
            # id-div step.  I think this is essentially "seqmagick convert
            # --include-from-file ..."
            sonar getfasta \
                -l ../../../../{input.seqids} \
                -f ../../../../{input.fasta} \
                -o ../../../../{output.fasta}
        """

def sonar_module_3_collect_inputs(wildcards):
    """List the inputs needed for SONAR module 3's collection step."""
    return igseq.sonar.sonar_islands_for_subject(
        SAMPLES, wildcards.subject, wildcards.chain, wildcards.chain_type)

def sonar_module_3_collect_seqs(wildcards, inputs):
    args = ["--seqs ../../../{key} --labels {val}".format() for key, val in inputs.items()]
    return " ".join(args)

rule sonar_module_3_collect:
    """Gather the selected island sequences from each specimen for a given subject and chain."""
    output:
        collected="analysis/sonar/{subject}/{chain}.{chain_type}/longitudinal/output/sequences/nucleotide/longitudinal-collected.fa"
    input: unpack(sonar_module_3_collect_inputs)
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd="analysis/sonar/{subject}/{chain}.{chain_type}/longitudinal",
        seqs=sonar_module_3_collect_seqs
    shell:
        """
            cd {params.wd}
            sonar merge_time {params.seqs}
        """

rule sonar_module_3_igphyml:
    """Run phylogenetic anaysis and generate tree across specimens for a given subject and chain."""
    output:
        tree="analysis/sonar/{subject}/{chain}.{chain_type}/longitudinal/output/longitudinal_igphyml.tree",
        inferred_nucl="analysis/sonar/{subject}/{chain}.{chain_type}/longitudinal/output/sequences/nucleotide/longitudinal_inferredAncestors.fa",
        inferred_prot="analysis/sonar/{subject}/{chain}.{chain_type}/longitudinal/output/sequences/amino_acid/longitudinal_inferredAncestors.fa",
        stats="analysis/sonar/{subject}/{chain}.{chain_type}/longitudinal/output/logs/longitudinal_igphyml_stats.txt"
    input:
        collected="analysis/sonar/{subject}/{chain}.{chain_type}/longitudinal/output/sequences/nucleotide/longitudinal-collected.fa",
        ref_hv="analysis/sonar/{subject}/{chain}.{chain_type}/germline.V.fasta",
        natives="analysis/sonar/{subject}/{chain}.{chain_type}/mab.fasta"
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        v_id=sonar_allele_v,
        args="-f"
    shell:
        """
            #set +e
            #set +u
            #source ~/miniconda3/bin/activate sonar
            cd $(dirname {input.collected})/../../..
            sonar igphyml \
                -v '{params.v_id}' \
                --lib ../germline.V.fasta \
                --natives {input.natives} \
                {params.args}
        """
