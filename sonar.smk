# A single SONAR project directory, as a format string.
# In general these are organized by specimen and equivalently timepoint but the
# longitudinal step gets its own project directory per-lineage.
WD_SONAR = Path("analysis/sonar/{subject}.{chain_type}/{specimen}")

import math
import csv
from collections import defaultdict

# if we have a custom version of the ID/DIV table, use that
ruleorder: sonar_module_2_id_div_island_alternate > sonar_module_2_id_div_island

# Make dynamically-defined rules for convenience, like "sonar_1_AJ09_gamma".
#
# The rules are:
#
#   * sonar_1_{subject}_{chain_type}
#   * sonar_2_{subject}_{chain_type}
#   * sonar_2_islands_{subject}_{chain_type}
#   * sonar_2_islands_{subject}_{chain_type}_{antibody_lineage}
#   * sonar_3_{subject}_{chain_type}
#   * sonar_3_{subject}_{chain_type}_{antibody_lineage}
#   * sonar_3_{subject}_{chain_type}_custom
#   * sonar_3_{subject}_{chain_type}_custom_{antibody_lineage}
#
# Get a full list of these rules and their descriptions with snakemake -l.
def sonar_setup_helper_rules():
    # specimen -> sets of chain_type (but only for IgG+)
    specimen_types = defaultdict(set)
    for attrs in SAMPLES.values():
        specimen_types[attrs["Specimen"]].add(attrs["Type"])
    # subject/chain_type -> sets of those specimens
    subject_type_map = defaultdict(set)
    for specimen, attrs in SPECIMENS.items():
        subject = attrs["Subject"]
        for chain_type in specimen_types[specimen]:
            if "IgG" in attrs["CellType"]:
                subject_type_map[(subject, chain_type)].add(specimen)
    # subject/chain_type -> sets of lineages
    subject_type_lineage_map = defaultdict(set)
    for key in subject_type_map:
        subject, chain_type = key
        for lineage, attrs in ANTIBODY_LINEAGES.items():
            if attrs["Subject"] == subject:
                subject_type_lineage_map[key].add(lineage)

    for key, specimens in subject_type_map.items():
        subject, chain_type = key
        # SONAR module 1 (annotation) for one subject and chain type
        targets = expand(
            "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_rearrangements.tsv",
            subject=subject, chain_type=chain_type, specimen=specimens)
        rule:
            f"SONAR Module 1 (annotation steps) for subject {subject} and chain type {chain_type}"
            name: f"sonar_1_{subject}_{chain_type}"
            input: targets
        # SONAR module 2 ID/DIV for one subject and chain type
        targets = expand(
            "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab",
            subject=subject, chain_type=chain_type, specimen=specimens)
        rule:
            f"SONAR Module 2 (ID/DIV table) for subject {subject} and chain type {chain_type}"
            name: f"sonar_2_{subject}_{chain_type}"
            input: targets
        # SONAR module 2 island for one subject and chain type and lineage
        # (rule for all lineages, and separate rules for each lineage)
        targets = expand("analysis/sonar/{subject}.{chain_type}/{specimen}/output/sequences/nucleotide/islandSeqs_{lineage}.fa",
            subject=subject, chain_type=chain_type, specimen=specimens,
            lineage=subject_type_lineage_map[key])
        rule:
            f"SONAR Module 2 (island selection) for subject {subject} and chain type {chain_type}, all lineages"
            name: f"sonar_2_islands_{subject}_{chain_type}"
            input: targets
        for lineage in subject_type_lineage_map[key]:
            targets = expand("analysis/sonar/{subject}.{chain_type}/{specimen}/output/sequences/nucleotide/islandSeqs_{lineage}.fa",
                subject=subject, chain_type=chain_type, specimen=specimens,
                lineage=lineage)
            rule:
                f"SONAR Module 2 (island selection) for subject {subject} and chain type {chain_type}, lineage {lineage}"
                name: f"sonar_2_islands_{subject}_{chain_type}_{lineage}"
                input: targets
    # SONAR module 3 (including tree PDF) for one subject and chain type
    for key, lineages in subject_type_lineage_map.items():
        if not lineages:
            continue
        # Automatic gathering and alignment across timepoints
        subject, chain_type = key
        targets = expand("analysis/sonar/{subject}.{chain_type}/longitudinal.auto.{lineage}/output/longitudinal.auto.{lineage}_igphyml.tree",
            subject=subject, chain_type=chain_type, lineage=lineages)
        rule:
            f"SONAR Module 3 (longitudinal analysis and tree) for subject {subject} and chain type {chain_type}, all lineages"
            name: f"sonar_3_{subject}_{chain_type}"
            input: targets
        for lineage in lineages:
            targets = expand("analysis/sonar/{subject}.{chain_type}/longitudinal.auto.{lineage}/output/longitudinal.auto.{lineage}_igphyml.tree",
                subject=subject, chain_type=chain_type, lineage=lineage)
            rule:
                f"SONAR Module 3 (longitudinal analysis and tree) for subject {subject} and chain type {chain_type}, lineage {lineage}"
                name: f"sonar_3_{subject}_{chain_type}_{lineage}"
                input: targets
        # Manually-prepared custom alignment input
        targets = expand("analysis/sonar/{subject}.{chain_type}/longitudinal.custom.{lineage}/output/longitudinal.custom.{lineage}_igphyml.tree",
            subject=subject, chain_type=chain_type, lineage=lineages)
        rule:
            f"SONAR Module 3 (longitudinal analysis and tree) for subject {subject} and chain type {chain_type}, all lineages, custom input"
            name: f"sonar_3_{subject}_{chain_type}_custom"
            input: targets
        for lineage in lineages:
            targets = expand("analysis/sonar/{subject}.{chain_type}/longitudinal.custom.{lineage}/output/longitudinal.custom.{lineage}_igphyml.tree",
                subject=subject, chain_type=chain_type, lineage=lineage)
            rule:
                f"SONAR Module 3 (longitudinal analysis and tree) for subject {subject} and chain type {chain_type}, lineage {lineage}, custom input"
                name: f"sonar_3_{subject}_{chain_type}_custom_{lineage}"
                input: targets

sonar_setup_helper_rules()

rule sonar_gather_mature:
    """Get heavy or light chain mature antibody sequences.

    This includes all lineages for the subject, if applicable, and all isolates
    that are not marked for exclusion from member identification.  This is
    stored in each project directory for a given subject so that you can
    optionally add custom sequences for a given timepoint.  (But, not at the
    top of the project directory, or SONAR will use it as input.)

    This will only keep the first occurrence of each sequence, so that mAb
    isolates that have identical sequences won't result in duplication of
    effort.
    """
    output: "analysis/sonar/{subject}.{chain_type}/{projdir}/mab/mab.fasta"
    run:
        if wildcards.chain_type in ["kappa", "lambda"]:
            seq_col = "LightSeq"
            light_locus = {"kappa": "IGK", "lambda": "IGL"}[wildcards.chain_type]
        else:
            seq_col = "HeavySeq"
        seen = {""} # always skip empty entries
        with open(output[0], "wt") as f_out:
            for seqid, attrs in ANTIBODY_ISOLATES.items():
                if seq_col == "LightSeq" and attrs["LineageAttrs"]["LightLocus"] != light_locus:
                    # Skip light chain sequences that are for the other locus
                    # than whatever was amplified here
                    continue
                if attrs.get("IncludeInTracing") != "Y":
                    # Skip any isolates marked for exclusion from the tracing
                    # itself
                    continue
                seq = attrs[seq_col]
                if attrs["LineageAttrs"]["Subject"] == wildcards.subject and seq not in seen:
                    f_out.write(f">{seqid} {attrs['Lineage']}\n")
                    f_out.write(attrs[seq_col]+"\n")
                    seen.add(seq)
        if seen == {""}:
            raise ValueError(
                f"No sequences found for subject {wildcards.subject} "
                f"and chain type {wildcards.chain_type}")

# I think this only comes up for including sequences for IgPhyML.  In this case
# we *will* keep duplicates so all mAb sequences will be in the tree, excluding
# those explicitly marked for exclusion for this step.
rule sonar_gather_mature_by_lineage:
    output: "analysis/sonar/{subject}.{chain_type}/{projdir}/mab/mab.{antibody_lineage}.fasta"
    run:
        if wildcards.chain_type in ["kappa", "lambda"]:
            seq_col = "LightSeq"
        else:
            seq_col = "HeavySeq"
        with open(output[0], "wt") as f_out:
            for seqid, attrs in ANTIBODY_ISOLATES.items():
                if attrs["Lineage"] == wildcards.antibody_lineage and attrs[seq_col] \
                        and attrs.get("IncludeInOutput") == "Y":
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

# Some of these SONAR tables and FASTAs are huge, so I'm compressing recent
# ones with xz, and adding this rule to default to compressed files if already
# available.  (Note to self: if you go and xz-compress some other module's
# stuff Snakemake will probably complain about ambiguous rules.)
ruleorder: sonar_module_1_decompress > sonar_module_1
rule sonar_module_1_decompress:
    output: (WD_SONAR/"output/{compressed}")
    input: (WD_SONAR/"output/{compressed}.xz")
    shell:
        """
            xz --decompress < {input} > {output}
        """

def input_sonar_germline(w):
    locus = {"kappa": "IGK", "lambda": "IGL"}.get(w.chain_type, "IGH")
    segs = ["V", "D", "J"] if locus == "IGH" else ["V", "J"]
    return {seg: f"analysis/germline/{w.subject}.{locus}/{seg}.fasta" for seg in segs}

rule sonar_module_1:
    output:
        fasta=WD_SONAR/"output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        rearr=WD_SONAR/"output/tables/{specimen}_rearrangements.tsv"
    input:
        # SONAR can use multiple input files and will automatically detect
        # them, but it requires them to be plaintext.
        unpack(input_sonar_germline),
        reads="analysis/samples-by-specimen/{specimen}.{chain_type}"
    log: WD_SONAR/"log.txt"
    singularity: "docker://jesse08/sonar"
    threads: 4
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        # 97% similarity was used in the SONAR vignette and Chaim says was more
        # appropriate back in the days of 454 sequencing, but they typically
        # use 99% now with the higher-quality Illumina sequencing.
        cluster_id_fract=config.get("sonar_cluster_id_fract", 0.99),
        cluster_min2=config.get("sonar_cluster_min2", 2),
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
        libd_arg=lambda w, input: "--dlib germline/D.fasta" if "D" in dict(input) else "--noD",
        input_v=lambda w, input: Path(input.V),
        input_d=lambda w, input: Path(input.D) if "D" in dict(input) else "",
        input_j=lambda w, input: Path(input.J),
        input_v_full=lambda w, input: Path(input.V).resolve(),
        input_d_full=lambda w, input: Path(input.D).resolve() if "D" in dict(input) else "",
        input_j_full=lambda w, input: Path(input.J).resolve(),
        jmotif=lambda w: JMOTIF[w.chain_type]
    shell:
        """
            for fqgz in {input.reads}/*.fastq.gz; do
                zcat $fqgz > {params.wd_sonar}/$(basename ${{fqgz%.gz}})
            done
            mkdir -p {params.wd_sonar}/germline
            cp {input.V} {input.J} {params.wd_sonar}/germline
            if [[ "{params.input_d}" != "" ]]; then
                cp {params.input_d} {params.wd_sonar}/germline
            fi
            (
              cd {params.wd_sonar}
              date
              echo "Running SONAR Module 1"
              echo
              echo "$(which sonar): $(sonar --version)"
              echo "Project directory: $PWD"
              echo "germline V: {params.input_v} -> {params.input_v_full}"
              if [[ "{params.input_d}" != "" ]]; then
                echo "germline D: {params.input_d} -> {params.input_d_full}"
              else
                echo "germline D: (No D)"
              fi
              echo "germline J: {params.input_j} -> {params.input_j_full}"
              echo "J motif: {params.jmotif}"
              echo "Cluster ID fract: {params.cluster_id_fract}"
              echo "Cluster min2: {params.cluster_min2}"
              echo
              set -x
              sonar blast_V --lib germline/V.fasta --derep --threads {threads} 2>&1
              sonar blast_J {params.libd_arg} --jlib germline/J.fasta --noC --threads {threads} 2>&1
              sonar finalize --jmotif '{params.jmotif}' --threads {threads} 2>&1
              sonar cluster_sequences --id {params.cluster_id_fract} --min2 {params.cluster_min2} 2>&1
            ) | tee -a {log}
        """

rule alternate_iddiv_with_igblast:
    # use IgBLAST results to supply the germline divergence values, rather than
    # SONAR's MUSCLE-derived values.
    output: iddiv=WD_SONAR/"output/tables/{specimen}_goodVJ_unique_id-div.alt.tab"
    input:
        iddiv=WD_SONAR/"output/tables/{specimen}_goodVJ_unique_id-div.tab",
        airr=lambda w: expand("analysis/igblast/custom-{{subject}}.{locus}/sonar/{{subject}}.{{chain_type}}/{{specimen}}/output/sequences/nucleotide/{{specimen}}_goodVJ_unique.fa.tsv.gz", locus={"kappa": "IGK", "lambda": "IGL"}.get(w.chain_type, "IGH"))
    run:
        with open(input.iddiv) as f_in, \
                gzip.open(str(input.airr), "rt") as f_in_airr, \
                open(output.iddiv, "wt") as f_out:
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
        unpack(input_sonar_germline),
        fasta=WD_SONAR / "output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        mab=WD_SONAR / "mab/mab.fasta"
    log: WD_SONAR/"log.id-div.txt"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_mab=lambda w, input: Path(input.mab).resolve(),
        # either "mismatch" (the default) or "ignore".  Note that this is used
        # for both identity and divergence calculations, and while ignoring
        # gaps can help with occasional mis-alignments with the V genes that
        # result in inflated divergence values, it also obscures key
        # differences in identity calculations (like for example with big gaps
        # in CDR3 alignments) that I think make gap=ignore infeasible to use.
        gap="mismatch"
    singularity: "docker://jesse08/sonar"
    threads: 4
    shell:
        """
            (
              cd {params.wd_sonar}
              date
              echo "Running SONAR Module 2 id-div"
              echo
              echo "$(which sonar): $(sonar --version)"
              echo "Project directory: $PWD"
              mab_seqs=$(grep -c "^>" {params.input_mab})
              echo "Reference antibodies: {input.mab} ($mab_seqs seqs)"
              echo "Gap setting: {params.gap}"
              mab={params.input_mab}
              set -x
              sonar id-div -g "germline/V.fasta" -a "$mab" -t {threads} --gap {params.gap}
            ) | tee -a {log}
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
        seen = {""} # same logic as for sonar_gather_mature
        with open(output[0], "wt") as f_out:
            for seqid, attrs in ANTIBODY_ISOLATES.items():
                seq = attrs[seq_col]
                if attrs["LineageAttrs"]["Subject"] == wildcards.subject \
                    and seq not in seen and \
                    attrs["Lineage"] == wildcards.antibody_lineage and \
                    attrs["IncludeInTracing"] == "Y":
                        f_out.write(f"{seqid}\n")
                        seen.add(seq)

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
    singularity: "docker://jesse08/sonar"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_iddiv=lambda w, input: Path(input.iddiv).resolve(),
        outprefix=lambda w, output: Path(output.seqids).stem
    shell:
        """
            mabargs=$(sed 's/^/--mab /' {input.mab})
            cd {params.wd_sonar}
            sonar get_island {params.input_iddiv} $mabargs --output {params.outprefix} --plotmethod binned
        """

rule sonar_module_2_id_div_island_alternate:
    output:
        seqids=WD_SONAR / "output/tables/islandSeqs_{antibody_lineage}.txt"
    input:
        iddiv=WD_SONAR / "output/tables/{specimen}_goodVJ_unique_id-div.alt.tab",
        mab=WD_SONAR / "mab/mab.{antibody_lineage}.txt"
    singularity: "docker://jesse08/sonar"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_iddiv=lambda w, input: Path(input.iddiv).resolve(),
        outprefix=lambda w, output: Path(output.seqids).stem
    shell:
        """
            mabargs=$(sed 's/^/--mab /' {input.mab})
            cd {params.wd_sonar}
            sonar get_island {params.input_iddiv} $mabargs --output {params.outprefix} --plotmethod binned
        """

rule sonar_module_2_id_div_getfasta:
    """SONAR 2: Extract FASTA matching selected island's seq IDs."""
    output:
        fasta=WD_SONAR / "output/sequences/nucleotide/islandSeqs_{txt}.fa"
    input:
        fasta=WD_SONAR / "output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        seqids=WD_SONAR / "output/tables/islandSeqs_{txt}.txt"
    singularity: "docker://jesse08/sonar"
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

# NOTE: specimens of the same timepoint will be labeled in the order they are
# given
def format_timepoints(specimens):
    # e.g.:
    # -12   "wkN12"
    #  -7   "wkN07"
    #   0   "wk000"
    #   4   "wk004"
    #   4   "wk004.2"
    #   8   "wk008"
    #  12   "wk012"
    try:
        # dict
        specs = specimens.values()
    except AttributeError:
        # list
        specs = specimens
    timepoints = [int(attrs["Timepoint"]) for attrs in specs]
    pads = [(tp<0) + math.floor(math.log10(max(1, abs(tp))) + 1) for tp in timepoints]
    padlen = max(pads)
    labels = ["wk" + str(num).replace("-", "N").zfill(padlen) for num in timepoints]
    trios = sorted(zip(timepoints, labels, specs), key=lambda trio: (trio[0], trio[1]))
    def tally_tp_label(label, specname, seen={}):
        if label in seen and specname not in seen[label]:
            seen[label].add(specname)
            out = label + "." + str(len(seen[label]))
            return out
        seen[label] = {specname}
        return label
    trios = [[tp, tally_tp_label(label, spec["Specimen"]), spec] for tp, label, spec in trios]
    return list(zip(*trios))

def sonar_module_3_collect_inputs(w):
    """List the inputs needed for SONAR module 3's collection step.

    This gives a dictionary of timepoint -> islandSeqs.fa pairs.
    """

    # first gather all relevant specimens.  This is the same logic as the
    # module 2 helper above, searching via samples to select specimens that
    # have amplicons for the right chain type.
    spec_names_here = set()
    for samp in SAMPLES.values():
        if samp["Type"] == w.chain_type and \
            samp["SpecimenAttrs"]["Subject"] == w.subject and \
            "IgG" in samp["SpecimenAttrs"]["CellType"]:
            spec_names_here.add(samp["Specimen"])
    spec_names_here = list(spec_names_here)
    if not spec_names_here:
        return []
    # sort specimens by timepoint and generate timepoints (integers) and
    # timepoint labels (strings)
    # first, an initial sort so any same-timepoint specimens are always labeled
    # the same.  Also specifically sort so that things like "...IGG" come
    # before "...IGG2"
    def spec_name_sorter(txt):
        match = re.search(r"IG[GM]([0-9]+)$", txt)
        if match:
            num = int(match.group(1))
            return (txt, num)
        return (txt, 0)
    spec_names_here.sort(key=spec_name_sorter)
    timepoints, labels, specs_here = format_timepoints([SPECIMENS[s] for s in spec_names_here])
    targets = expand(
        "analysis/sonar/{subject}.{chain_type}/{other_specimen}/"
        "output/sequences/nucleotide/islandSeqs_{antibody_lineage}.fa",
        subject=w.subject,
        chain_type=w.chain_type,
        other_specimen=[attrs["Specimen"] for attrs in specs_here],
        antibody_lineage=w.antibody_lineage)
    targets = {label: target for label, target in zip(labels, targets)}
    return targets

def sonar_module_3_collect_param_seqs(_, input):
    args = [" --labels {key} --seqs {val}".format(key=key, val=Path(val).resolve()) for key, val in input.items()]
    return " ".join(args)

rule sonar_module_3_collect:
    output:
        collected="analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}/output/sequences/nucleotide/longitudinal.{word}.{antibody_lineage}-collected.fa"
    input: unpack(sonar_module_3_collect_inputs)
    singularity: "docker://jesse08/sonar"
    threads: 4
    params:
        wd_sonar=lambda w: expand("analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}", **w),
        seqs=sonar_module_3_collect_param_seqs
    shell:
        """
            cd {params.wd_sonar}
            sonar merge_time {params.seqs}
        """

def sonar_module_3_igphyml_param_v_id(wildcards):
    if dict(wildcards).get("word") == "custom":
        return ""
    key_default = "germline_v"
    key = "germline_v_{subject}_{chain_type}_{antibody_lineage}".format(**wildcards)
    v_call = config.get(key, config.get(key_default))
    if v_call:
        return v_call
    if wildcards.chain_type in ["kappa", "lambda"]:
        v_call = ANTIBODY_LINEAGES[wildcards.antibody_lineage]["VL"]
    else:
        v_call = ANTIBODY_LINEAGES[wildcards.antibody_lineage]["VH"]
    if not v_call:
        raise ValueError("No V assigned for %s" % wildcards)
    return v_call

def input_for_sonar_module_3_igphyml(w):
    proj = "longitudinal.{word}.{antibody_lineage}"
    sonar_dir = "analysis/sonar/{subject}.{chain_type}/" + proj + "/"
    if w.word == "auto":
        targets = input_sonar_germline(w)
        targets.update({
            "collected": sonar_dir + "output/sequences/nucleotide/" + proj + "-collected.fa",
            "natives": sonar_dir + "mab/mab.{antibody_lineage}.fasta"
        })
    else:
        targets = {
            "alignment": "analysis/sonar/{subject}.{chain_type}/alignment.{antibody_lineage}.{word}.fa"
        }
    return targets

# (Setting this up as a checkpoint means we can have downstream rules structured
# according to the contents of these output files, for example, what ancestor
# sequences are inferred.)
checkpoint sonar_module_3_igphyml:
    """SONAR 3: Run IgPhyML with auto or custom alignment and generate tree across specimens.

    For the automatic alignment case ("auto" as the keyword before the lineage
    name), this will build a tree with the appropriate germline V sequence for
    the lineage included and will root the tree on that V sequence.

    For the custom alignment case (any other keyword), this will given the
    first sequence ID in the alignment as the --root for sonar igphyml.
    """
    output:
        tree="analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}/output/longitudinal.{word}.{antibody_lineage}_igphyml.tree",
        inferred_nucl="analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}/output/sequences/nucleotide/longitudinal.{word}.{antibody_lineage}_inferredAncestors.fa",
        inferred_prot="analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}/output/sequences/amino_acid/longitudinal.{word}.{antibody_lineage}_inferredAncestors.fa",
        stats="analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}/output/logs/longitudinal.{word}.{antibody_lineage}_igphyml_stats.txt",
        # alignment produced here if word=auto, just copied from input (so the output is consistent) otherwise
        alignment="analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}/work/phylo/longitudinal.{word}.{antibody_lineage}_aligned.afa"
    input: unpack(input_for_sonar_module_3_igphyml)
    log: "analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}/log.txt"
    singularity: "docker://jesse08/sonar"
    params:
        wd_sonar=lambda w: expand("analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}", **w),
        seed=config.get("sonar_igphyml_seed", 123),
        aln_type=lambda w: "auto" if w.word == "auto" else "custom",
        args="-f",
        input_germline_v=lambda w, input: Path(input.V).resolve() if w.word == "auto" else "",
        input_natives=lambda w, input: Path(input.natives).resolve() if w.word == "auto" else "",
        input_alignment=lambda w, input: Path(input.alignment).resolve() if w.word != "auto" else "",
        output_alignment=lambda w, output: Path(output.alignment).resolve(),
        v_id=sonar_module_3_igphyml_param_v_id,
    shell:
        """
            # See sonar_module_3_igphyml_auto about LANG
            unset LANG
            (
              cd {params.wd_sonar}
              date
              echo "Running SONAR Module 3 igphyml"
              echo
              echo "Alignment type: {params.aln_type}"
              echo "$(which sonar): $(sonar --version)"
              echo "Project directory: $PWD"
              echo "Random seed: {params.seed}"
              if [[ {params.aln_type} == "auto" ]]; then
                echo "Germline V ID for root: {params.v_id}"
                echo "Germline V library: {params.input_germline_v}"
                echo "Natives FASTA: {params.input_natives}"
                set -x
                sonar igphyml \
                    -v '{params.v_id}' \
                    --lib {params.input_germline_v} \
                    --natives {params.input_natives} \
                    --seed {params.seed} \
                    {params.args} 2>&1
              else
                root=$(head -n 1 {params.input_alignment} | cut -c 2- | cut -f 1 -d ' ')
                echo "Detected seq ID for tree root from first FASTA record: $root"
                set -x
                sonar igphyml --root "$root" -i {params.input_alignment} --seed {params.seed} {params.args} 2>&1
                cp {params.input_alignment} {params.output_alignment}
              fi

            ) | tee -a {log}
        """

rule sonar_make_natives_table:
    output:
        table="analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}/natives.tab"
    params:
        # we don't actually need the files here but this has the logic to name
        # the timepoints
        tp_fastas=sonar_module_3_collect_inputs
    log: "analysis/sonar/{subject}.{chain_type}/longitudinal.{word}.{antibody_lineage}/log.sonar_make_natives_table.txt"
    run:
        keys = list(dict(params.tp_fastas).keys())
        mabs = {}
        if wildcards.chain_type in ["kappa", "lambda"]:
            seq_col = "LightSeq"
        else:
            seq_col = "HeavySeq"
        with open(log[0], "w") as log_out:
            log_out.write(f"starting timepoints {keys}\n")
            log_out.write(f"lineage {wildcards.antibody_lineage}\n")
            log_out.write(f"seq_col {seq_col}\n")
            # TODO don't assume that the isolates come later than all the
            # bulk-sequenced timepoints.  That'll require a bit of reorganizing
            # to how this is done, though.
            for seqid, attrs in ANTIBODY_ISOLATES.items():
                if attrs["Lineage"] == wildcards.antibody_lineage and attrs[seq_col]:
                    tp = "wk" + attrs["Timepoint"]
                    if tp not in mabs:
                        mabs[tp] = []
                    if tp not in keys:
                        keys.append(tp)
                    mabs[tp].append(seqid)
                    log_out.write(f"{tp} -> {seqid}\n")
            with open(output.table, "wt") as f_out:
                for key in keys:
                    if key in mabs:
                        for mab in mabs[key]:
                            f_out.write(f"{key}\t{mab}\t{mab}\n")
                    else:
                        f_out.write(f"{key}\t\t\n")

rule sonar_module_3_draw_tree:
    output:
        # here "thing" can include the lineage field since the lineage isn't
        # directly used at this point
        tree_img="analysis/sonar/{subject}.{chain_type}/longitudinal{thing}/output/longitudinal{thing}_igphyml.tree{suffix,.*}.pdf"
    input:
        tree="analysis/sonar/{subject}.{chain_type}/longitudinal{thing}/output/longitudinal{thing}_igphyml.tree",
        natives_tab="analysis/sonar/{subject}.{chain_type}/longitudinal{thing}/natives{suffix}.tab"
    singularity: "docker://jesse08/sonar"
    params:
        wd_sonar=lambda w: expand("analysis/sonar/{subject}.{chain_type}/longitudinal{thing}", **w),
        # (Jumping through some hoops to run inside of the SONAR project
        # directory to avoid leaving things like SONAR_command_history.log
        # lying around at the top level)
        input_tree=lambda w, input: Path(input.tree).resolve(),
        input_natives_tab=lambda w, input: Path(input.natives_tab).resolve(),
        output_tree_img=lambda w, input, output: Path(output.tree_img).resolve(),
    # Running via xvfb-run since the ETE toolkit requires X11 to render for
    # some reason
    shell:
        """
            cd {params.wd_sonar}
            xvfb-run sonar display_tree \
                -t {params.input_tree} \
                -n {params.input_natives_tab} \
                --showAll \
                -o {params.output_tree_img}
        """
