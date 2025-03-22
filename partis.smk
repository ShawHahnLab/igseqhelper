"""
Run partis on all IgG+ material across specimens for one subject and chain type
(using SONAR's clustered sequences), and partition all of those sequences plus
the corresponding mature isolate sequences.
"""

def partis_specimens(subject, chain_type):
    """List specimen names that should be used for partis for one subject+chain type"""
    specs = []
    for spec_attrs in SPECIMENS.values():
        if spec_attrs["Subject"] == subject and "IgG+" in spec_attrs["CellType"]:
            specs.append(spec_attrs["Specimen"])
    for samp_attrs in SAMPLES.values():
        spec_attrs = samp_attrs.get("SpecimenAttrs", {})
        if samp_attrs["Type"] == chain_type and \
                spec_attrs.get("Subject") == subject and \
                "IgG+" in spec_attrs.get("CellType"):
            specs.append(spec_attrs["Specimen"])
    return specs

def input_for_partis_ngs_fasta(w):
    """input paths for partis_ngs_fasta, as dict with timepoint label as key"""
    spec_names_here = partis_specimens(w.subject, w.chain_type)
    timepoints, labels, specs_here = format_timepoints([SPECIMENS[s] for s in spec_names_here])
    targets = expand(
        "analysis/sonar/{subject}.{chain_type}/{specimen}/output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        subject = w.subject,
        chain_type = w.chain_type,
        specimen = [attrs["Specimen"] for attrs in specs_here])
    targets = {label: target for label, target in zip(labels, targets)}
    targets["excludes"] = f"analysis/partis/{w.subject}.{w.chain_type}/ngs.excludes.fasta"
    return targets

rule partis_ngs_fasta_excludes:
    """Placeholder for excludes from NGS sequences given to partis"""
    # partis crashes sometimes with really weird input sequences.  Use this to
    # exclude them.
    output: "analysis/partis/{subject}.{chain_type}/ngs.excludes.fasta"
    shell: "touch {output}"

rule partis_ngs_fasta:
    """Cross-timepoint NGS seqs for partis"""
    output: "analysis/partis/{subject}.{chain_type}/ngs.fasta"
    input: unpack(input_for_partis_ngs_fasta)
    run:
        seq_col = "LightSeq" if wildcards.chain_type in ["kappa", "lambda"] else "HeavySeq"
        excludes = {rec.id: str(rec.seq) for rec in SeqIO.parse(input.excludes, "fasta")}
        with open(output[0], "w") as f_out:
            for label, fasta in input.items():
                # Tried to have a separate input for the excludes, but
                # Snakemake says "Cannot combine named input file (name
                # timepoints) with unpack()", so evidently I need a workaround:
                # just skip excludes as it's handled above, and these should
                # all be the per-timepoint SONAR FASTA inputs
                if label == "excludes":
                    continue
                for rec in SeqIO.parse(fasta, "fasta"):
                    rec.description = ""
                    rec.id = label + "-" + rec.id
                    if rec.id in excludes:
                        if excludes[rec.id] != str(rec.seq):
                            raise ValueError(f"Seq ID/content mismatch from excludes: {rec.id}")
                    else:
                        SeqIO.write(rec, f_out, "fasta-2line")

rule partis_isolates_fasta:
    """Cross-lineage isolate seqs for partis"""
    output: "analysis/partis/{subject}.{chain_type}/isolates.fasta"
    run:
        seq_col = "LightSeq" if wildcards.chain_type in ["kappa", "lambda"] else "HeavySeq"
        with open(output[0], "w") as f_out:
            for isolate_attrs in ANTIBODY_ISOLATES.values():
                if isolate_attrs[seq_col]:
                    lineage_attrs = isolate_attrs.get("LineageAttrs", {})
                    isolate_prefix = re.sub("-.*", "", isolate_attrs["Isolate"])
                    subject = lineage_attrs.get("Subject")
                    # include an isolate for this subject if the metadata
                    # specifically links the two, or failing that, if the
                    # isolate name begins with the subject name
                    if (subject and subject == wildcards.subject) or \
                            (not subject and isolate_prefix == wildcards.subject):
                        seqid = isolate_attrs["Isolate"]
                        seq = isolate_attrs[seq_col]
                        f_out.write(f">{seqid}\n{seq}\n")

rule partis_combo_fasta:
    """Combined isolate+NGS seqs for partis"""
    output: "analysis/partis/{subject}.{chain_type}/combined.fasta"
    input:
        isolates="analysis/partis/{subject}.{chain_type}/isolates.fasta",
        ngs="analysis/partis/{subject}.{chain_type}/ngs.fasta"
    shell: "cat {input} > {output}"

rule partis_cache_params:
    """Run partis cache-parameters on a cross-timepoint NGS FASTA"""
    output: directory("analysis/partis/{subject}.{chain_type}/params")
    input: "analysis/partis/{subject}.{chain_type}/ngs.fasta"
    log:
        main="analysis/partis/{subject}.{chain_type}/cache_params.log.txt",
        conda="analysis/partis/{subject}.{chain_type}/cache_params.conda_build.txt"
    params:
        locus=lambda w: {"kappa": "igk", "lambda": "igl"}.get(w.chain_type, "igh"),
        species=config.get("partis_species", "macaque"),
        seed=config.get("partis_random_seed", 1),
        # I have the dependencies provided via conda but the actual partis
        # install still lives in a big ball of stuff in one directory,
        # unfortunately.  The first like of the shell commands will ensure this
        # variable is set.
        partis=os.getenv("PARTIS_HOME", "")
    conda: "envs/partis.yaml"
    threads: 14
    shell:
        """
            if [[ ! -n "{params.partis}" ]]; then
              echo "Need path to partis install"
              exit 1
            fi
            (
              date
              echo "PARTIS_HOME: {params.partis}"
              echo "locus: {params.locus}"
              echo "species: {params.species}"
              echo "random seed: {params.seed}"
            ) >> {log.main}
            {params.partis}/bin/partis cache-parameters --n-procs {threads} \
              --random-seed {params.seed} --locus {params.locus} --species {params.species} \
              --infname {input} --parameter-dir {output} 2>&1 | tee -a {log.main}
        """

rule partis_partition:
    """Run partis partition on a combined isolate+NGS FASTA with existing cached parameters"""
    output: "analysis/partis/{subject}.{chain_type}/partitions.yaml"
    input:
        params_dir="analysis/partis/{subject}.{chain_type}/params",
        fasta="analysis/partis/{subject}.{chain_type}/combined.fasta"
    log:
        main="analysis/partis/{subject}.{chain_type}/partition.log.txt",
        conda="analysis/partis/{subject}.{chain_type}/partition.conda_build.txt"
    params:
        locus=lambda w: {"kappa": "igk", "lambda": "igl"}.get(w.chain_type, "igh"),
        species=config.get("partis_species", "macaque"),
        seed=config.get("partis_random_seed", 1),
        partis=os.getenv("PARTIS_HOME", "")
    conda: "envs/partis.yaml"
    threads: 14
    shell:
        """
            if [[ ! -n "{params.partis}" ]]; then
              echo "Need path to partis install"
              exit 1
            fi
            conda list --explicit > {log.conda}
            (
              date
              echo "PARTIS_HOME: {params.partis}"
              echo "locus: {params.locus}"
              echo "species: {params.species}"
              echo "random seed: {params.seed}"
            ) >> {log.main}
            {params.partis}/bin/partis partition --n-procs {threads} --no-naive-vsearch \
              --random-seed {params.seed} --locus {params.locus} --species {params.species} \
              --infname {input.fasta} --parameter-dir {input.params_dir} --outfname {output} 2>&1 | tee -a {log.main}
        """

rule partis_partition_airr:
    """Convert Partis' YAML output for partition data to AIRR TSV"""
    output: "analysis/partis/{subject}.{chain_type}/partitions.airr.tsv"
    input: "analysis/partis/{subject}.{chain_type}/partitions.yaml"
    log:
        main="analysis/partis/{subject}.{chain_type}/partition.log.txt",
        conda="analysis/partis/{subject}.{chain_type}/partition.conda_build.txt"
    params:
        partis=os.getenv("PARTIS_HOME", "")
    conda: "envs/partis.yaml"
    shell:
        """
            if [[ ! -n "{params.partis}" ]]; then
              echo "Need path to partis install"
              exit 1
            fi
            conda list --explicit > {log.conda}
            (
              date
              echo "PARTIS_HOME: {params.partis}"
            ) >> {log.main}
            {params.partis}/bin/parse-output.py --airr-output {input} {output} 2>&1 | tee -a {log.main}
        """

### Below is really reporting logic, but, no time for proper organization here

def input_for_partis_seq_lineage_info(w):
    if w.chain_type in ("kappa", "lambda"):
        raise ValueError("light chain not supported")
    # always required the final AIRR table as input
    path = Path(expand("analysis/partis/{subject}.{chain_type}", **w)[0])
    # ...and IgBLAST's AIRR too
    igblast_path = Path(expand("analysis/igblast/custom-{subject}.IGH/partis/{subject}.{chain_type}", **w)[0])
    targets = {
        "airr": path/"partitions.airr.tsv",
        "airr_igblast": igblast_path/"combined.fasta.tsv.gz",
        "isolates": ancient("metadata/isolates.csv")}
    # If there's a CSV provided with NGS sequence lineage info, use that also,
    # so we can assign custom lineage IDs to NGS seqs
    path_annotations = path/"ngs_lineages.csv"
    if path_annotations.exists():
        targets["ngs_annots"] = path_annotations
    return targets

rule partis_seq_lineage_info:
    """Report sequence and lineage info from partis, merging metadata for our seqs+isolates"""
    output: "analysis/partis/{subject}.{chain_type}/seq_lineages.csv"
    input: unpack(input_for_partis_seq_lineage_info)
    run:
        cmd = "partis_seq_lineage_info.py {input.airr} {output} -i {input.isolates} -A {input.airr_igblast} --all"
        if "ngs_annots" in dict(input):
            cmd += " -n {input.ngs_annots}"
        shell(cmd)

rule partis_lineages:
    """Summarize partis info per-lineage-group, one row per group+timepoint+category"""
    output: "analysis/partis/{subject}.{chain_type}/lineage_groups.csv"
    input: "analysis/partis/{subject}.{chain_type}/seq_lineages.csv"
    shell: "partis_lineages.py {input} {output}"

rule partis_lineages_summary:
    """Summarize partis info per-lineage-group further, one row per lineage group"""
    output: "analysis/partis/{subject}.{chain_type}/lineage_groups_summary.csv"
    input: "analysis/partis/{subject}.{chain_type}/lineage_groups.csv"
    shell: "partis_lineages_summary.py {input} {output}"
