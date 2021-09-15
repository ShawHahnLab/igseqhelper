"""
SONAR processing across the three modules.

There's one interactive step here, the island selection from the
identity/divergence plots.  This requires X11 to show the interactive plot.
Most rules here require --use-singularity in snakemake.

This whole setup is pretty brittle and would probably take some work to run
successfully anywhere else.
"""

import igseq.sonar
from igseq.data import MetadataError, transpose_sample_md

SAMPLE_MD_IGG = transpose_sample_md(SAMPLES, "IgG+")

# The working directory (as a format string) for a single SONAR analysis.  One
# level up will contain information shared by specimens but unique to that
# lineage.  One more, shared by lineages but unique for that
# subject/chain/chain type.
WD_SONAR = Path("analysis/sonar/{subject}/{chain}.{chain_type}/{antibody_lineage}/{specimen}")

def setup_sonar_targets(sample_md_igg, antibody_lineages):
    """Create dictionary of targets for various SONAR rules.

    The interleaved paths with many templated variables get really intricate
    here, so I'm hiding the complexity in this setup fuction and storing the
    targets in a single unified dictionary based on SONAR module.
    """
    sonar_combos = igseq.sonar.setup_sonar_combos(sample_md_igg, antibody_lineages)
    pattern_root = "analysis/sonar/{subject}/{chain}.{chain_type}/{antibody_lineage}"
    patterns = {
        "prep": pattern_root + "/{specimen}/{specimen}.fastq",
        "module_1": pattern_root + "/{specimen}/output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        "module_2_auto": pattern_root + "/{specimen}/output/tables/{specimen}_lineages.txt",
        "module_2_id_div": pattern_root + "/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab",
        "module_2_id_div_island": pattern_root + "/{specimen}/output/sequences/nucleotide/{specimen}_islandSeqs.fa",
        "module_3": pattern_root + "/longitudinal/output/longitudinal_igphyml.tree"}
    targets = {key: expand(pattern, zip, **sonar_combos) for key, pattern in patterns.items()}
    return targets

TARGETS_SONAR = setup_sonar_targets(SAMPLE_MD_IGG, ANTIBODY_LINEAGES)

rule all_sonar_prep_input:
    input: TARGETS_SONAR["prep"]

rule all_sonar_module_1:
    input: TARGETS_SONAR["module_1"]

rule all_sonar_module_2_auto:
    input: TARGETS_SONAR["module_2_auto"]

rule all_sonar_module_2_id_div:
    input: TARGETS_SONAR["module_2_id_div"]

rule all_sonar_module_2_id_div_island:
    input: TARGETS_SONAR["module_2_id_div_island"]

rule all_sonar_module_3:
    input: TARGETS_SONAR["module_3"]

def sonar_gather_germline_inputs(wildcards):
    """Get IgDiscover per-specimen paths for a given per-subject output."""
    # IgG+ will have gamma for heavy, but the germline will have come from mu.
    return igseq.sonar.igdiscover_final_db(
        SAMPLES, wildcards.subject, wildcards.chain, wildcards.chain_type, wildcards.segment)

rule sonar_prep_input_from_presto:
    """Gather input sequences from pRESTO's initial processing."""
    output: (WD_SONAR / "{specimen}.fastq")
    input: "analysis/presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq"
    params:
        input=lambda w, input: Path(input[0]).resolve()
    shell: "ln {params.input} {output}"

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
    output: WD_SONAR.parent.parent / "germline.{segment}.fasta"
    input: sonar_gather_germline_inputs
    params:
        output=lambda w, input, output: Path(output[0]).resolve()
    run:
        igseq.sonar.gather_germline(input, params.output)

rule sonar_module_1:
    """SONAR 1: Collect good sequences and tabulate AIRR info from raw input."""
    output:
        fasta=(WD_SONAR / "output/sequences/nucleotide/{specimen}_goodVJ_unique.fa"),
        rearr=(WD_SONAR / "output/tables/{specimen}_rearrangements.tsv")
    input: unpack(igseq.sonar.sonar_module_1_inputs)
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_fastq=lambda w, input: Path(input.fastq).resolve(),
        # 97% similarity was used in the SONAR vignette and Chaim says was more
        # appropriate back in the days of 454 sequencing, but they typically
        # use 99% now with the higher-quality Illumina sequencing.
        cluster_id_fract=.97,
        cluster_min2=2,
        # What we give for the V(D)J command-line argments and the J motif
        # depends on which chain we're doing.  We could handle this during the
        # rule itself with a run: directive but then it won't let us use
        # singularity.
        # Note that if we DON'T give an explicit argument for --jmotif, SONAR
        # will use a regular expression combining the default heavy AND light
        # chain patterns.  Most of the time this is fine, but in some cases (so
        # far I saw this only in IGHJ6*01) the gene just happens to have a
        # short section matching the wrong chain's motif, and the CDR3 region
        # is marked incorrectly.  Instead we'll always explicitly specify
        # either the heavy or light motif to avoid this problem.
        # (SONAR says --jmotif specifies "Conserved nucleotide sequence
        # indicating the start of FWR4 on the J gene.  Defaults to either TGGGG
        # for heavy chains or TT[C|T][G|A]G for light chains; set manually for
        # custom light chain libraries or for species with a different motif.")
        libv_arg=lambda w, input: "--lib " + str(Path(input.V).resolve()),
        libd_arg=lambda w, input: "D" in input and "--dlib " + str(Path(input.D).resolve()) or "--noD",
        libj_arg=lambda w, input: "--jlib " + str(Path(input.J).resolve()),
        jmotif=lambda w: {"heavy": "TGGGG", "light": "TT[C|T][G|A]G"}[w.chain]
    shell:
        """
            cd {params.wd_sonar}
            sonar blast_V --fasta {params.input_fastq} {params.libv_arg} --derep --threads {threads}
            sonar blast_J {params.libd_arg} {params.libj_arg} --noC --threads {threads}
            sonar finalize --noclean --jmotif '{params.jmotif}' --threads {threads}
            sonar cluster_sequences --id {params.cluster_id_fract} --min2 {params.cluster_min2}
        """

rule sonar_gather_mature:
    """Get heavy or light chain mature antibody sequences for a lineage."""
    output: WD_SONAR.parent / "mab.fasta"
    run:
        igseq.sonar.gather_mature(
            ANTIBODY_ISOLATES,
            wildcards.antibody_lineage,
            wildcards.chain, 
            output[0])

### Non-interactive aspect of module 2

def sonar_allele_v(wildcards):
    """Get the sequence ID for the expected germline V allele"""
    seqid = igseq.sonar.get_antibody_allele(
        ANTIBODY_LINEAGES, wildcards.antibody_lineage, wildcards.subject, wildcards.chain, "V")
    return seqid

def sonar_allele_j(wildcards):
    """Get the sequence ID for the expected germline J allele"""
    seqid = igseq.sonar.get_antibody_allele(
        ANTIBODY_LINEAGES, wildcards.antibody_lineage, wildcards.subject, wildcards.chain, "J")
    return seqid

# Automated intradonor analysis and grouping
rule sonar_module_2:
    """SONAR 2: Automated intradonor analysis.

    We're not using this currently; see id_div rules below.
    """ 
    output: str(WD_SONAR / "output/tables/{specimen}_lineages.txt")
    input: unpack(igseq.sonar.sonar_module_2_inputs)
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        v_id=sonar_allele_v,
        j_id=sonar_allele_j
    shell:
        """
            # Start from the top-level directory for this specimen's analysis
            cd {params.wd_sonar}
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
rule sonar_module_2_id_div:
    """SONAR 2: Calc seqs' identity to MABs vs divergence from germline V."""
    output:
        iddiv=WD_SONAR / "output/tables/{specimen}_goodVJ_unique_id-div.tab"
    input:
        fasta=WD_SONAR / "output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        mab=WD_SONAR.parent / "mab.fasta",
        germline_v=WD_SONAR.parent.parent / "germline.V.fasta"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_v=lambda w, input: Path(input.germline_v).resolve(),
        input_mab=lambda w, input: Path(input.mab).resolve(),
        # either "mismatch" (the default) or "ignore".  Note that this is used
        # for both identity and divergence calculations, and while ignoring
        # gaps can help with occasional mis-alignments with the V genes that
        # result in inflated divergence values, it also obscures key
        # differences in identity calculations (like for example with big gaps
        # in CDR3 alignments) that I think make gap=ignore infeasible to use.
        gap="mismatch"
    singularity: "docker://scharch/sonar"
    threads: 4
    shell:
        """
            cd {params.wd_sonar}
            libv={params.input_v}
            mab={params.input_mab}
            sonar id-div -g "$libv" -a "$mab" -t {threads} --gap {params.gap}
        """

# Part 2 of 3: Select specific clusters from id/div plots
# NOTE this step is interactive over X11
rule sonar_module_2_id_div_island:
    """SONAR 2: Select "island" of interest from ID/DIV plots."""
    output:
        seqids=WD_SONAR / "output/tables/islandSeqs.txt"
    input:
        iddiv=WD_SONAR / "output/tables/{specimen}_goodVJ_unique_id-div.tab"
    singularity: "docker://scharch/sonar"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_iddiv=lambda w, input: Path(input.iddiv).resolve(),
        mab="=a", # =a means use all antibodies in mab input file
    shell:
        """
            cd {params.wd_sonar}
            sonar get_island {params.input_iddiv} --mab "{params.mab}"
        """

# Part 3 of 3: Extract those goodVJ cluster sequences given in the id/div lists
rule sonar_module_2_id_div_getfasta:
    """SONAR 2: Extract FASTA matching selected island's seq IDs."""
    output:
        fasta=WD_SONAR / "output/sequences/nucleotide/{specimen}_islandSeqs.fa"
    input:
        fasta=WD_SONAR / "output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        seqids=WD_SONAR / "output/tables/islandSeqs.txt"
    singularity: "docker://scharch/sonar"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_fasta=lambda w, input: Path(input.fasta).resolve(),
        input_seqids=lambda w, input: Path(input.seqids).resolve(),
        output_fasta=lambda w, input, output: Path(output.fasta).resolve()
    shell:
        """
            cd {params.wd_sonar}
            # Pull out sequences using our selected island from the above
            # id-div step.  I think this is essentially "seqmagick convert
            # --include-from-file ..."
            sonar getFastaFromList \
                -l {params.input_seqids} \
                -f {params.input_fasta} \
                -o {params.output_fasta}
        """

def sonar_module_3_collect_inputs(wildcards):
    """List the inputs needed for SONAR module 3's collection step.

    This gives a dictionary of timepoint -> islandSeqs.fa pairs.  We need to
    include all other specimens that include samples for the same chain type as
    specified here, so this searches for samples and then works backwards to
    specimens.
    """

    pattern = str(WD_SONAR.parent /
        "{other_specimen}/output/sequences/nucleotide/{other_specimen}_islandSeqs.fa")
    def filter_sample_md(vec):
        match = lambda trio: trio[1] == wildcards.subject and trio[2] == wildcards.chain_type
        trios = zip(vec, SAMPLE_MD_IGG["subjects"], SAMPLE_MD_IGG["chaintypes"])
        return [trio[0] for trio in trios if match(trio)]
    timepoints = filter_sample_md(SAMPLE_MD_IGG["timepoints"])
    specimens = filter_sample_md(SAMPLE_MD_IGG["specimens"])
    targets = expand(pattern,
        subject=wildcards.subject,
        chain=wildcards.chain,
        chain_type=wildcards.chain_type,
        antibody_lineage=wildcards.antibody_lineage,
        other_specimen=specimens)
    # zero-pad the timepoints, so for example we end up with ["08", "12", ...
    # instead of ["12", "8", ...
    padlen = max([len(txt) for txt in timepoints])
    timepoints = [txt.zfill(padlen) for txt in timepoints]
    # What I've done so far will give repeated entries, so we'll get just
    # unique ones now.
    pairs = list(set(zip(timepoints, targets)))
    pairs = sorted(pairs)
    # Allow for repeated timepoints (for the rare case of multiple specimens
    # per timepoint) by appending .2, .3, etc.
    # (note to self: be careful with this mutable default argument.)
    def tally_tp(tp, seen={}):
        if tp in seen:
            seen[tp] += 1
            out = tp + "." + str(seen[tp])
            return out
        seen[tp] = 1
        return tp
    targets = {"wk" + tally_tp(tp): target for tp, target in pairs}
    return targets

def sonar_module_3_collect_seqs(wildcards, input):
    args = [" --labels {key} --seqs {val}".format(key=key, val=Path(val).resolve()) for key, val in input.items()]
    return " ".join(args)

rule sonar_module_3_collect:
    """SONAR 3: Gather selected sequences from each specimen for given subject and chain."""
    output:
        collected=WD_SONAR / "output/sequences/nucleotide/longitudinal-collected.fa"
    input: unpack(sonar_module_3_collect_inputs)
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        seqs=sonar_module_3_collect_seqs
    shell:
        """
            cd {params.wd_sonar}
            sonar merge_time {params.seqs}
        """

rule sonar_module_3_igphyml:
    """SONAR 3: Run phylogenetic analysis and generate tree across specimens."""
    output:
        tree=WD_SONAR / "output/longitudinal_igphyml.tree",
        inferred_nucl=WD_SONAR / "output/sequences/nucleotide/longitudinal_inferredAncestors.fa",
        inferred_prot=WD_SONAR / "output/sequences/amino_acid/longitudinal_inferredAncestors.fa",
        stats=WD_SONAR / "output/logs/longitudinal_igphyml_stats.txt"
    input:
        collected=WD_SONAR / "output/sequences/nucleotide/longitudinal-collected.fa",
        germline_v=WD_SONAR.parent.parent / "germline.V.fasta",
        natives=WD_SONAR.parent / "mab.fasta"
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_germline_v=lambda w, input: Path(input.germline_v).resolve(),
        input_natives=lambda w, input: Path(input.natives).resolve(),
        v_id=sonar_allele_v,
        args="-f"
    shell:
        """
            # Singularity passes through environment variables by default,
            # which in my setup includes LANG=en_US.UTF-8 which is unrecognized
            # by the container, which then causes a few warnings to be written
            # to stderr by a perl script called by SONAR, which then triggers
            # SONAR to crash when it sees the non-empty stderr.
            # So yeah let's just unset the LANG.
            unset LANG
            cd {params.wd_sonar}
            sonar igphyml \
                -v '{params.v_id}' \
                --lib {params.input_germline_v} \
                --natives {params.input_natives} \
                {params.args}
        """

rule sonar_module_3_draw_tree:
    """SONAR 3: Render a tree visualization from the IgPhyML results.

    natives.tab is a tab-separated table of metadata described in SONAR's
    display_tree.py.  Currently this has to be supplied manually.
    """
    output:
        tree_img=WD_SONAR / "output/longitudinal_igphyml.tree.pdf"
    input:
        tree=WD_SONAR / "output/longitudinal_igphyml.tree",
        natives_tab=WD_SONAR / "natives.tab"
    singularity: "docker://scharch/sonar"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_tree=lambda w, input: Path(input.tree).resolve(),
        input_natives_tab=lambda w, input: Path(input.natives_tab).resolve(),
        output_tree_img=lambda w, input, output: Path(output.tree_img).resolve()
    shell:
        """
            sonar display_tree \
                -t {params.input_tree} \
                -n {params.input_natives_tab} \
                --showAll \
                -o {params.output_tree_img}
        """
