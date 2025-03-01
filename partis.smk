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
            [ -n "{params.partis}" ]
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
            [ -n "{params.partis}" ]
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
            [ -n "{params.partis}" ]
            conda list --explicit > {log.conda}
            (
              date
              echo "PARTIS_HOME: {params.partis}"
            ) >> {log.main}
            {params.partis}/bin/parse-output.py --airr-output {input} {output} 2>&1 | tee -a {log.main}
        """

###

def input_for_partis_seq_lineage_info(w):
    # always required the final AIRR table as input
    path = Path(expand("analysis/partis/{subject}.{chain_type}", **w)[0])
    targets = {
        "airr": path/"partitions.airr.tsv",
        "isolates": ancient("metadata/isolates.csv")}
    # If there's a CSV provided with NGS sequence lineage info, use that also,
    # so we can assign custom lineage IDs to NGS seqs
    path_annotations = path/"ngs_lineages.csv"
    if path_annotations.exists():
        targets["ngs_annots"] = path_annotations
    return targets

rule partis_seq_lineage_info:
    """Make CSV of partis and my lineage assignment info for all lineages including our isolates"""
    # This should probably be a reporting rule, really, but, starting here for
    # now
    output: "analysis/partis/{subject}.{chain_type}/seq_lineages.csv"
    input: unpack(input_for_partis_seq_lineage_info)
    run:
        if "ngs_annots" in dict(input):
            shell("partis_seq_lineage_info.py {input.airr} {output} -i {input.isolates} -n {input.ngs_annots}")
        else:
            shell("partis_seq_lineage_info.py {input.airr} {output} -i {input.isolates}")

rule partis_lineages:
    # summarize per-lineage-group
    output: "analysis/partis/{subject}.{chain_type}/lineage_groups.csv"
    input: "analysis/partis/{subject}.{chain_type}/seq_lineages.csv"
    run:
        with open(input[0]) as f_in:
            seq_info = list(DictReader(f_in))
        groups = defaultdict(list)
        for row in seq_info:
            groups[row["lineage_group"]].append(row)
        timepoints = set()
        categories = set()
        keys_out = set()
        out = []
        for lineage_group, rows in groups.items():
            # tally seqs per category and timepoint
            totals = defaultdict(int)
            for row in rows:
                timepoint = int(row["timepoint"])
                category = row["category"]
                timepoints.add(timepoint)
                categories.add(category)
                key = (category, timepoint)
                totals[key] += 1
            # flatten into one dictionary per lineage group
            row_out = {"lineage_group": lineage_group}
            per_category_totals = defaultdict(int)
            for category, timepoint in totals:
                key_out = f"{category}_wk{timepoint}"
                keys_out.add(key_out)
                per_category_totals[category] += totals[(category, timepoint)]
                row_out[key_out] = totals[(category, timepoint)]
            for category, num in per_category_totals.items():
                row_out[f"{category}_total"] = num
            out.append(row_out)
        # sort output rows by decreasing totals
        def sorter(row):
            total = sum(num or 0 for key, num in row.items() if key.endswith("_total"))
            return (-total, row["lineage_group"])
        out.sort(key=sorter)
        # figure out output columns and their order
        cols = ["lineage_group"]
        for timepoint in sorted(timepoints):
            if timepoint != "total":
                for category in sorted(categories):
                    key_out = f"{category}_wk{timepoint}"
                    if key_out in keys_out:
                        cols.append(key_out)
        for category in sorted(categories):
            cols.append(f"{category}_total")
        with open(output[0], "w") as f_out:
            writer = DictWriter(f_out, cols, lineterminator="\n")
            writer.writeheader()
            writer.writerows(out)
