# The working directory (as a format string) for a single SONAR analysis.
WD_SONAR = Path("analysis/sonar/{subject}.{chain_type}/{specimen}")

def input_helper_sonar(w, pattern):
    # Take all specimens for this subject and the corresponding amplicons.
    # IgG+ is implicit in these rules but other types can be requested
    # manually.
    specimens = set()
    for samp in SAMPLES.values():
        if samp["Type"] == w.chain_type and \
            samp["SpecimenAttrs"]["Subject"] == w.subject and \
            "IgG" in samp["SpecimenAttrs"]["CellType"]:
            specimens.add(samp["Specimen"])
    return expand(
        pattern,
        subject=w.subject,
        chain_type=w.chain_type,
        specimen=specimens)

# For one subject and amplicon, process all specimens for module 1
rule helper_sonar_module_1_by_subject:
    output: touch("analysis/sonar/{subject}.{chain_type}/module1.done")
    input: lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_rearrangements.tsv")

# For one lineage and amplicon, process all specimens for module 2 (just ID/DIV)
rule helper_sonar_module_2_by_lineage:
    output: touch("analysis/sonar/{subject}.{chain_type}/module2.done")
    input: lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab")

# For one lineage and amplicon, make all the island FASTAs via the plots and lists
rule helper_sonar_module_2_island_by_lineage:
    output: touch("analysis/sonar/{subject}.{chain_type}/module2island.done")
    input: lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/sequences/nucleotide/{specimen}_islandSeqs.fa")

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
    output: WD_SONAR/"mab/mab.fasta"
    run:
        if wildcards.chain_type in ["kappa", "lambda"]:
            seq_col = "LightSeq"
        else:
            seq_col = "HeavySeq"
        with open(output[0], "wt") as f_out:
            for seqid, attrs in ANTIBODY_ISOLATES.items():
                if attrs["AntibodyLineageAttrs"]["Subject"] == wildcards.subject and attrs[seq_col]:
                    f_out.write(f">{seqid}\n")
                    f_out.write(attrs[seq_col])

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
rule sonar_module_2_id_div_island:
    output:
        seqids=WD_SONAR / "output/tables/islandSeqs_{antibody_lineage}.txt"
    input:
        iddiv=WD_SONAR / "output/tables/{specimen}_goodVJ_unique_id-div.tab",
        mab=WD_SONAR / "mab/mab.{antibody_lineage}.txt"
    singularity: "docker://scharch/sonar"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_iddiv=lambda w, input: Path(input.iddiv).resolve()
    shell:
        """
            mabargs=$(sed 's/^/--mab /' {input.mbargs})
            cd {params.wd_sonar}
            sonar get_island {params.input_iddiv} $mbargs
        """
