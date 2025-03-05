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
    """Make CSV of partis and my lineage assignment info for all lineages including our isolates"""
    # This should probably be a reporting rule, really, but, starting here for
    # now
    output: "analysis/partis/{subject}.{chain_type}/seq_lineages.csv"
    input: unpack(input_for_partis_seq_lineage_info)
    run:
        cmd = "partis_seq_lineage_info.py {input.airr} {output} -i {input.isolates} -A {input.airr_igblast} --all"
        if "ngs_annots" in dict(input):
            cmd += " -n {input.ngs_annots}"
        shell(cmd)

rule partis_lineages:
    # summarize per-lineage-group, one row per group+timepoint+category
    output: "analysis/partis/{subject}.{chain_type}/lineage_groups.csv"
    input: "analysis/partis/{subject}.{chain_type}/seq_lineages.csv"
    run:
        cols = [
            "lineage_group",
            "v_family",
            "v_identity_min",
            "v_identity_max",
            "d_call",
            "junction_aa_v_min",
            "junction_aa_v_max",
            "junction_aa_length_min",
            "junction_aa_length_max",
            "timepoint",
            "category",
            "total"]
        with open(input[0]) as f_in:
            seq_info = list(DictReader(f_in))
        groups = defaultdict(list)
        for row in seq_info:
            # group rows by lineage group, including those with none assigned
            # (in case all sequences were included in the input CSV)
            groups[row["lineage_group"]].append(row)
        out = []
        for lineage_group, rows in groups.items():
            # within each lineage group, group rows by timepoint and category
            chunk = defaultdict(list)
            for row in rows:
                timepoint = int(row["timepoint"])
                category = row["category"]
                chunk[(timepoint, category)].append(row)
            for key, rows in chunk.items():
                v_family = ""
                d_call = ""
                junction_aa_length_min = min(int(row["junction_aa_length"]) for row in rows)
                junction_aa_length_max = max(int(row["junction_aa_length"]) for row in rows)
                juncts = None
                if lineage_group:
                    v_family = "/".join(sorted({row["v_family"] for row in rows}))
                    d_call = set()
                    for row in rows:
                        d_call = d_call | set(re.split("[/,]", row["d_call"]))
                    d_call = "/".join(sorted(d_call - {""}))
                    junction_aa_length_min = min(int(row["junction_aa_length"]) for row in rows)
                    junction_aa_length_max = max(int(row["junction_aa_length"]) for row in rows)
                    juncts = [(float(row["v_identity"]) if row["v_identity"] else 0, row["junction_aa"]) for row in rows]
                    juncts.sort()
                # set up output per group+timepoint+category
                out.append({
                    "lineage_group": lineage_group,
                    "v_family": v_family,
                    "v_identity_min": f"{juncts[0][0]:.2f}" if juncts else "",
                    "v_identity_max": f"{juncts[-1][0]:.2f}" if juncts else "",
                    "d_call": d_call,
                    "junction_aa_v_min": juncts[0][1] if juncts else "",
                    "junction_aa_v_max": juncts[-1][1] if juncts else "",
                    "junction_aa_length_min": junction_aa_length_min,
                    "junction_aa_length_max": junction_aa_length_max,
                    "timepoint": key[0],
                    "category": key[1],
                    "total": len(rows)})
        # total counts per lineage goup across timepoints and categories
        totals = defaultdict(int)
        for row in out:
            totals[row["lineage_group"]] += row["total"]
        # sort by lineage group by decreasing total count, then by increasing
        # timepoint, then by category
        def sorter(row):
            return (
                -totals[row["lineage_group"]],
                row["lineage_group"],
                row["timepoint"],
                row["category"])
        out.sort(key=sorter)
        with open(output[0], "w") as f_out:
            writer = DictWriter(f_out, cols, lineterminator="\n")
            writer.writeheader()
            writer.writerows(out)

rule partis_lineages_summary:
    # summarize per-lineage-group, wide format with one row per group
    output: "analysis/partis/{subject}.{chain_type}/lineage_groups_summary.csv"
    input: "analysis/partis/{subject}.{chain_type}/lineage_groups.csv"
    run:
        with open(input[0]) as f_in:
            info = list(DictReader(f_in))
        # define categories, timepoints, output columns
        categories = sorted({row["category"] for row in info})
        timepoints = sorted({int(row["timepoint"]) for row in info})
        cols = ["lineage_group"]
        for category in categories:
            for timepoint in timepoints:
                cols.append(f"total_{category}_wk{timepoint}")
            cols.append(f"total_{category}")
        cols += [
            "v_family",
            "v_identity_min",
            "v_identity_max",
            "d_call",
            "junction_aa_v_min",
            "junction_aa_v_max",
            "junction_aa_length_min",
            "junction_aa_length_max"]
        # group by lineage group
        groups = defaultdict(list)
        for row in info:
            groups[row["lineage_group"]].append(row)
        # define output rows
        out = []
        for group, rows in groups.items():
            # D call(s) across everything (parsing and re-formatting any with
            # more than one to properly handle ["X/Y", "X"] -> "X/Y").  If
            # there are more than four, throw up our hands and just put "???"
            d_call = set()
            for row in rows:
                d_call = d_call | set(row["d_call"].split("/"))
            d_call = d_call - {""}
            if len(d_call) > 4:
                d_call = "???"
            else:
                d_call = "/".join(sorted(d_call))
            # V family across everything (should only be one!  but if not, same
            # idea as for D, just no upper limit.)
            v_family = set()
            for row in rows:
                v_family = v_family | set(row["v_family"].split("/"))
            v_family = v_family - {""}
            v_family = "/".join(sorted(v_family))
            # Minimum and maximum junction AA length across all sequences
            # across all categories and timepoints
            junction_aa_length_min = min(int(row["junction_aa_length_min"]) for row in rows)
            junction_aa_length_max = max(int(row["junction_aa_length_max"]) for row in rows)
            # The pair of junction AA sequences associated with the lowest V
            # identity and the highest V identity across all sequences across
            # all categories and timepoints
            juncts_min = [(float(row["v_identity_min"]) if row["v_identity_min"] else 0, row["junction_aa_v_min"]) for row in rows]
            juncts_min.sort()
            juncts_max= [(float(row["v_identity_max"]) if row["v_identity_max"] else 0, row["junction_aa_v_max"]) for row in rows]
            juncts_max.sort(reverse=True)
            row_out = {
                "lineage_group": group,
                "v_family": v_family,
                "v_identity_min": juncts_min[0][0],
                "v_identity_max": juncts_max[0][0],
                "d_call": d_call,
                "junction_aa_v_min": juncts_min[0][1],
                "junction_aa_v_max": juncts_max[0][1],
                "junction_aa_length_min": junction_aa_length_min,
                "junction_aa_length_max": junction_aa_length_max,
                }
            totals = defaultdict(int)
            for row in rows:
                category = row["category"]
                timepoint = row["timepoint"]
                row_out[f"total_{category}_wk{timepoint}"] = int(row["total"])
                totals[category] += int(row["total"])
            for category in categories:
                row_out[f"total_{category}"] = totals[category]
            out.append(row_out)
        with open(output[0], "w") as f_out:
            writer = DictWriter(f_out, cols, lineterminator="\n")
            writer.writeheader()
            writer.writerows(out)
