# A single SONAR project directory, as a format string.
# In general these are organized by specimen and equivalently timepoint but the
# longitudinal step gets its own project directory per-lineage.
WD_SONAR = Path("analysis/sonar/{subject}.{chain_type}/{specimen}")
WD_SONAR_LONG = Path("analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}")

import csv

# if we have a custom version of the ID/DIV table, use that
ruleorder: sonar_module_2_id_div_island_alternate > sonar_module_2_id_div_island

def input_helper_sonar(w, pattern):
    # Take all specimens for this subject and the corresponding amplicons.
    # IgG+ is implicit in these rules but other types can be requested
    # manually.
    parts = vars(w)
    specimens = set()
    subject = parts.get("subject")
    if not subject and "antibody_lineage" in parts:
        subject = ANTIBODY_LINEAGES[w.antibody_lineage]["Subject"]
    if not subject:
        raise ValueError
    for samp in SAMPLES.values():
        if samp["Type"] == w.chain_type and \
            "IgG" in samp["SpecimenAttrs"]["CellType"]:
            if samp["SpecimenAttrs"]["Subject"] == subject:
                specimens.add(samp["Specimen"])
    parts["subject"] = subject
    parts["specimen"] = specimens
    return expand(pattern, **parts)

# For one subject and amplicon, process all specimens for module 1
rule helper_sonar_module_1_by_subject:
    output: touch("analysis/sonar/{subject}.{chain_type}/module1.done")
    input: lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_rearrangements.tsv")

# For one amplicon, process all specimens for module 2 (just ID/DIV, for members of all lineages)
rule helper_sonar_module_2:
    output: touch("analysis/sonar/{subject}.{chain_type}/module2.done")
    input: lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab")

# For one lineage and amplicon, make all the island FASTAs via the plots and lists
rule helper_sonar_module_2_island_by_lineage:
    output: touch("analysis/sonar/{subject}.{chain_type}/module2island.{antibody_lineage}.done")
    input: lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/sequences/nucleotide/islandSeqs_{antibody_lineage}.fa")

# For one lienage and amplicon, make the IgPhyML tree and related outputs
rule helper_sonar_module_3_igphyml:
    output: touch("analysis/sonar/{subject}.{chain_type}/module3tree.{antibody_lineage}.done")
    input: lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/output/longitudinal-{antibody_lineage}_igphyml.tree")

# just a default setup using SONARRamesh -> IgDiscover results
# D will be ignored for light chian but is always there (which makes the
# snakemake rules easy)
# for light chain we want the corresponding light chain from the IgDiscover
# files, but for heavy chain it'd be mu.
IGG_IGM = {
    "alpha": "mu",
    "delta": "mu",
    "gamma": "mu",
    "mu": "mu",
    "epsilon": "mu",
    "kappa": "kappa",
    "lambda": "lambda"}
rule sonar_germline:
    output: WD_SONAR.parent/"germline.{segment}.fasta"
    input: lambda w: expand("analysis/igdiscover/SONARRamesh/{chain_type}/{{subject}}/final/database/{{segment}}.fasta", chain_type=IGG_IGM[w.chain_type])
    shell: "cp {input} {output}"

rule sonar_gather_mature:
    """Get heavy or light chain mature antibody sequences.

    This includes all lineages for the subject, if applicable.  This is stored
    in each project directory for a given subject so that you can optionally
    add custom sequences for a given timepoint.  (But, not at the top of the
    proejct directory, or SONAR will use it as input.)
    """
    output: "analysis/sonar/{subject}.{chain_type}/{projdir}/mab/mab.fasta"
    run:
        if wildcards.chain_type in ["kappa", "lambda"]:
            seq_col = "LightSeq"
        else:
            seq_col = "HeavySeq"
        with open(output[0], "wt") as f_out:
            for seqid, attrs in ANTIBODY_ISOLATES.items():
                if attrs["AntibodyLineageAttrs"]["Subject"] == wildcards.subject and attrs[seq_col]:
                    f_out.write(f">{seqid}\n")
                    f_out.write(attrs[seq_col]+"\n")

rule sonar_gather_mature_by_lineage:
    output: "analysis/sonar/{subject}.{chain_type}/{projdir}/mab/mab.{antibody_lineage}.fasta"
    run:
        if wildcards.chain_type in ["kappa", "lambda"]:
            seq_col = "LightSeq"
        else:
            seq_col = "HeavySeq"
        with open(output[0], "wt") as f_out:
            for seqid, attrs in ANTIBODY_ISOLATES.items():
                if attrs["AntibodyLineage"] == wildcards.antibody_lineage and attrs[seq_col]:
                    f_out.write(f">{seqid}\n")
                    f_out.write(attrs[seq_col]+"\n")

# Any heavy chain will use a simple constant j motif, but light chains have a
# variety of possible sequences.
JMOTIF = {
    "alpha": "TGGGG",
    "delta": "TGGGG",
    "gamma": "TGGGG",
    "mu": "TGGGG",
    "epsilon": "TGGGG",
    "kappa": "TT[C|T][G|A]G",
    "lambda": "TT[C|T][G|A]G"}

rule sonar_module_1:
    output:
        fasta=(WD_SONAR/"output/sequences/nucleotide/{specimen}_goodVJ_unique.fa"),
        rearr=(WD_SONAR/"output/tables/{specimen}_rearrangements.tsv")
    input:
        # SONAR can use multiple input files and will automatically detect
        # them, but it requires them to be plaintext.
        reads="analysis/samples-by-specimen/{specimen}.{chain_type}",
        # Using the already-prepped VDJ reference files from the sonar_germline
        # rule above
        V="analysis/sonar/{subject}.{chain_type}/germline.V.fasta",
        D="analysis/sonar/{subject}.{chain_type}/germline.D.fasta",
        J="analysis/sonar/{subject}.{chain_type}/germline.J.fasta"
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        # 97% similarity was used in the SONAR vignette and Chaim says was more
        # appropriate back in the days of 454 sequencing, but they typically
        # use 99% now with the higher-quality Illumina sequencing.
        cluster_id_fract=.99,
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
        jmotif=lambda w: JMOTIF[w.chain_type]
    shell:
        """
            for fqgz in {input.reads}/*.fastq.gz; do
                zcat $fqgz > {params.wd_sonar}/$(basename ${{fqgz%.gz}})
            done
            cd {params.wd_sonar}
            sonar blast_V {params.libv_arg} --derep --threads {threads}
            sonar blast_J {params.libd_arg} {params.libj_arg} --noC --threads {threads}
            sonar finalize --noclean --jmotif '{params.jmotif}' --threads {threads}
            sonar cluster_sequences --id {params.cluster_id_fract} --min2 {params.cluster_min2}
        """

rule igblast_sonar:
    # custom igblast with one of SONAR's FASTA files, giving AIRR-format TSV
    # output.
    output: WD_SONAR/"output/tables/{thing}.igblast.tsv"
    input:
        fasta=WD_SONAR/"output/sequences/nucleotide/{thing}.fa",
        germline=expand("{prefix}/germline.{segment}.fasta", prefix=WD_SONAR.parent, segment=["V", "D", "J"])
    threads: 8
    shell:
        """
            igseq igblast -t {threads} -S rhesus -r {input.germline} -Q {input.fasta} -outfmt 19 -out {output}
        """

rule alternate_iddiv_with_igblast:
    # use IgBLAST results to supply the germline divergence values, rather than
    # SONAR's MUSCLE-derived values.
    output: iddiv=WD_SONAR/"output/tables/{specimen}_goodVJ_unique_id-div.alt.tab"
    input:
        iddiv=WD_SONAR/"output/tables/{specimen}_goodVJ_unique_id-div.tab",
        airr=WD_SONAR/"output/tables/{specimen}_goodVJ_unique.igblast.tsv"
    run:
        pass
        with open(input.iddiv) as f_in, open(input.airr) as f_in_airr, open(output.iddiv, "wt") as f_out:
            iddiv = csv.DictReader(f_in, delimiter="\t")
            airr = csv.DictReader(f_in_airr, delimiter="\t")
            writer = csv.DictWriter(f_out, fieldnames=iddiv.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for iddiv_row, airr_row in zip(iddiv, airr):
                iddiv_id = iddiv_row["sequence_id"]
                airr_id = airr_row["sequence_id"]
                if iddiv_id != airr_id:
                    raise ValueError(
                        "row mismatch between ID-DIV table and AIRR: %s vs %s", iddiv_id, airr_id)
                sonar_gene = iddiv_row["v_gene"]
                sonar_div = float(iddiv_row["germ_div"])
                airr_gene = airr_row["v_call"].split("*")[0]
                airr_div = 100 - float(airr_row["v_identity"])
                iddiv_row["germ_div"] = f"{airr_div:.2f}"
                writer.writerow(iddiv_row)

### Manually-guided aspect of Module 2
# Identity/Divergence plots and island identification
# The first part-- calculating the identity/divergence coordinates-- is
# automatic.  The second part-- selecting points on the scatterplots-- is
# manual and requires X11. (The third is just extracting the selected sequences
# from the original FASTA.)
rule sonar_module_2_id_div:
    output:
        iddiv=WD_SONAR / "output/tables/{specimen}_goodVJ_unique_id-div.tab"
    input:
        fasta=WD_SONAR / "output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        mab=WD_SONAR / "mab/mab.fasta",
        germline_v=WD_SONAR.parent / "germline.V.fasta"
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

rule sonar_list_members_for_lineage:
    # gets all the seq IDs as in the input mab antibody FASTA and save in a
    # text file, so they can be given to sonar get_island like "--mab ID1 --mab
    # ID2 ..." (as opposed to "--mab =a" for all of them from the ID/DIV table)
    output: WD_SONAR / "mab/mab.{antibody_lineage}.txt"
    run:
        if wildcards.chain_type in ["kappa", "lambda"]:
            seq_col = "LightSeq"
        else:
            seq_col = "HeavySeq"
        with open(output[0], "wt") as f_out:
            for seqid, attrs in ANTIBODY_ISOLATES.items():
                if attrs["AntibodyLineageAttrs"]["Subject"] == wildcards.subject \
                    and attrs[seq_col] and \
                    attrs["AntibodyLineage"] == wildcards.antibody_lineage:
                        f_out.write(f"{seqid}\n")

# NOTE this step is interactive over X11
##
# To get X working with Snakemake + singularity I had to append both of these
# arguments to the snakemake call:
#
#     --use-singularity --singularity-args "-B /home -B /data/home -H /home/$USER"
#
# I think these home directory options are needed so that the Xauthority stuff
# works but I'm not totally sure.   Using a regular Singularity image doesn't
# need this so I think it must be something about how Snakemake calls
# singularity.
rule sonar_module_2_id_div_island:
    output:
        seqids=WD_SONAR / "output/tables/islandSeqs_{antibody_lineage}.txt"
    input:
        iddiv=WD_SONAR / "output/tables/{specimen}_goodVJ_unique_id-div.tab",
        mab=WD_SONAR / "mab/mab.{antibody_lineage}.txt"
    singularity: "docker://scharch/sonar"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_iddiv=lambda w, input: Path(input.iddiv).resolve(),
        outprefix=lambda w, output: Path(output.seqids).stem
    shell:
        """
            mabargs=$(sed 's/^/--mab /' {input.mab})
            cd {params.wd_sonar}
            sonar get_island {params.input_iddiv} $mabargs --output {params.outprefix}
        """

rule sonar_module_2_id_div_island_alternate:
    output:
        seqids=WD_SONAR / "output/tables/islandSeqs_{antibody_lineage}.txt"
    input:
        iddiv=WD_SONAR / "output/tables/{specimen}_goodVJ_unique_id-div.alt.tab",
        mab=WD_SONAR / "mab/mab.{antibody_lineage}.txt"
    singularity: "docker://scharch/sonar"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_iddiv=lambda w, input: Path(input.iddiv).resolve(),
        outprefix=lambda w, output: Path(output.seqids).stem
    shell:
        """
            mabargs=$(sed 's/^/--mab /' {input.mab})
            cd {params.wd_sonar}
            sonar get_island {params.input_iddiv} $mabargs --output {params.outprefix}
        """

rule sonar_module_2_id_div_getfasta:
    """SONAR 2: Extract FASTA matching selected island's seq IDs."""
    output:
        fasta=WD_SONAR / "output/sequences/nucleotide/islandSeqs_{txt}.fa"
    input:
        fasta=WD_SONAR / "output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        seqids=WD_SONAR / "output/tables/islandSeqs_{txt}.txt"
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

def sonar_module_3_collect_inputs(w):
    """List the inputs needed for SONAR module 3's collection step.

    This gives a dictionary of timepoint -> islandSeqs.fa pairs.
    """

    # first gather all relevant specimens.  This is the same logic as the
    # module 2 helper above, searching via samples to select specimens that
    # have amplicons for the right chain type.
    specimens = set()
    for samp in SAMPLES.values():
        if samp["Type"] == w.chain_type and \
            samp["SpecimenAttrs"]["Subject"] == w.subject and \
            "IgG" in samp["SpecimenAttrs"]["CellType"]:
            specimens.add(samp["Specimen"])
    specimens = list(specimens)
    # infer input FASTAs from those specimens
    targets = expand(
        "analysis/sonar/{subject}.{chain_type}/{other_specimen}/"
        "output/sequences/nucleotide/islandSeqs_{antibody_lineage}.fa",
        subject=w.subject,
        chain_type=w.chain_type,
        other_specimen=specimens,
        antibody_lineage=w.antibody_lineage)
    # get timepoints and zero-pad, so for example we end up with ["08", "12",
    # ...  instead of ["12", "8", ...
    timepoints = [SPECIMENS[spec]["Timepoint"] for spec in specimens]
    padlen = max([len(txt) for txt in timepoints])
    timepoints = [txt.zfill(padlen) for txt in timepoints]
    # Pair up timepoints and target FASTAs and sort by timepoint
    pairs = sorted(zip(timepoints, targets))
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

def sonar_module_3_collect_param_seqs(_, input):
    args = [" --labels {key} --seqs {val}".format(key=key, val=Path(val).resolve()) for key, val in input.items()]
    return " ".join(args)

# analysis/sonar/T681.gamma/longitudinal-T681/output/sequences/nucleotide/longitudinal-T681-collected.fa
rule sonar_module_3_collect:
    output:
        collected=WD_SONAR_LONG / "output/sequences/nucleotide/longitudinal-{antibody_lineage}-collected.fa"
    input: unpack(sonar_module_3_collect_inputs)
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR_LONG), **w),
        seqs=sonar_module_3_collect_param_seqs
    shell:
        """
            cd {params.wd_sonar}
            sonar merge_time {params.seqs}
        """

def sonar_module_3_igphyml_param_v_id(wildcards):
    if wildcards.chain_type in ["kappa", "lambda"]:
        return ANTIBODY_LINEAGES[wildcards.antibody_lineage]["VL"]
    return ANTIBODY_LINEAGES[wildcards.antibody_lineage]["VH"]

rule sonar_module_3_igphyml:
    """SONAR 3: Run phylogenetic analysis and generate tree across specimens."""
    output:
        tree=WD_SONAR_LONG / "output/longitudinal-{antibody_lineage}_igphyml.tree",
        inferred_nucl=WD_SONAR_LONG / "output/sequences/nucleotide/longitudinal-{antibody_lineage}_inferredAncestors.fa",
        inferred_prot=WD_SONAR_LONG / "output/sequences/amino_acid/longitudinal-{antibody_lineage}_inferredAncestors.fa",
        stats=WD_SONAR_LONG / "output/logs/longitudinal-{antibody_lineage}_igphyml_stats.txt"
    input:
        collected=WD_SONAR_LONG / "output/sequences/nucleotide/longitudinal-{antibody_lineage}-collected.fa",
        germline_v=WD_SONAR_LONG.parent / "germline.V.fasta",
        natives=WD_SONAR_LONG / "mab/mab.{antibody_lineage}.fasta"
    singularity: "docker://scharch/sonar"
    threads: 4
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR_LONG), **w),
        input_germline_v=lambda w, input: Path(input.germline_v).resolve(),
        input_natives=lambda w, input: Path(input.natives).resolve(),
        v_id=sonar_module_3_igphyml_param_v_id,
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

rule sonar_make_natives_table:
    output:
        table=WD_SONAR_LONG/"natives.tab"
    # we don't actually need the files here but this has the logic to name the
    # timepoints
    input: unpack(sonar_module_3_collect_inputs)
    run:
        keys = list(dict(input).keys())
        mabs = {}
        if wildcards.chain_type in ["kappa", "lambda"]:
            seq_col = "LightSeq"
        else:
            seq_col = "HeavySeq"
        for seqid, attrs in ANTIBODY_ISOLATES.items():
            if attrs["AntibodyLineage"] == wildcards.antibody_lineage and attrs[seq_col]:
                tp = "wk" + attrs["Timepoint"]
                if tp not in mabs:
                    mabs[tp] = []
                if tp not in keys:
                    keys.append(tp)
                mabs[tp].append(seqid)
        with open(output.table, "wt") as f_out:
            for key in keys:
                if key in mabs:
                    for mab in mabs[key]:
                        f_out.write(f"{key}\t{mab}\t{mab}\n")
                else:
                    f_out.write(f"{key}\t\t\n")

rule sonar_module_3_draw_tree:
    output:
        tree_img=WD_SONAR_LONG / "output/longitudinal-{antibody_lineage}_igphyml.tree.pdf"
    input:
        tree=WD_SONAR_LONG / "output/longitudinal-{antibody_lineage}_igphyml.tree",
        natives_tab=WD_SONAR_LONG / "natives.tab"
    singularity: "docker://scharch/sonar"
    # Running via xvfb-run since the ETE toolkit requires X11 to render for
    # some reason
    shell:
        """
            xvfb-run sonar display_tree \
                -t {input.tree} \
                -n {input.natives_tab} \
                --showAll \
                -o {output.tree_img}
        """
