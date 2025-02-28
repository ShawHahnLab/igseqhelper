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
    return targets

rule partis_ngs_fasta:
    """Cross-timepoint NGS seqs for partis"""
    output: "analysis/partis/{subject}.{chain_type}/ngs.fasta"
    input: unpack(input_for_partis_ngs_fasta)
    run:
        seq_col = "LightSeq" if wildcards.chain_type in ["kappa", "lambda"] else "HeavySeq"
        with open(output[0], "w") as f_out:
            for label, fasta in input.items():
                for rec in SeqIO.parse(fasta, "fasta"):
                    rec.description = ""
                    rec.id = label + "-" + rec.id
                    SeqIO.write(rec, f_out, "fasta-2line")

rule partis_isolates_fasta:
    """Cross-lineage isolate seqs for partis"""
    output: "analysis/partis/{subject}.{chain_type}/isolates.fasta"
    run:
        seq_col = "LightSeq" if wildcards.chain_type in ["kappa", "lambda"] else "HeavySeq"
        with open(output[0], "w") as f_out:
            for isolate_attrs in ANTIBODY_ISOLATES.values():
                if isolate_attrs["LineageAttrs"]["Subject"] == wildcards.subject and isolate_attrs[seq_col]:
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
        locus="igh",
        species="macaque",
        # I have the dependencies provided via conda but the actual partis
        # install still lives in a big ball of stuff in one directory,
        # unfortunately.  The first like of the shell commands will ensure this
        # variable is set.
        partis=os.getenv("PARTIS_HOME", "")
    conda: "envs/partis.yaml"
    threads: 14
    shell:
        """
            [ -n "{params.partis}" ]
            (
              date
              echo "PARTIS_HOME: {params.partis}"
              echo "locus: {params.locus}"
              echo "species: {params.species}"
            ) >> {log.main}
            {params.partis}/bin/partis cache-parameters --n-procs {threads} \
              --locus {params.locus} --species {params.species} \
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
        locus="igh",
        species="macaque",
        partis=os.getenv("PARTIS_HOME", "")
    conda: "envs/partis.yaml"
    threads: 14
    shell:
        """
            [ -n "{params.partis}" ]
            conda list --explicit > {log.conda}
            (
              date
              echo "PARTIS_HOME: {params.partis}"
              echo "locus: {params.locus}"
              echo "species: {params.species}"
            ) >> {log.main}
            {params.partis}/bin/partis partition --n-procs {threads} \
              --no-naive-vsearch --locus {params.locus} --species {params.species} \
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
            [ -n "{params.partis}" ]
            conda list --explicit > {log.conda}
            (
              date
              echo "PARTIS_HOME: {params.partis}"
            ) >> {log.main}
            {params.partis}/bin/parse-output.py --airr-output {input} {output} 2>&1 | tee -a {log.main}
        """
