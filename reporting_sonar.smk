import itertools
from collections import defaultdict
from statistics import median

# dynamic rules for per-lineage targets
for antibody_lineage, attrs in ANTIBODY_LINEAGES.items():
    chain_types = ["gamma"]
    try:
        chain_types.append({"L": "lambda", "K": "kappa"}[attrs["VL"][2]])
    except:
        pass
    things = [
        "island_stats_summary.csv",
        "rearrangements_summary.tsv",
        "members.csv",
        "igphyml_collected.csv",
        "igphyml_ancestors.csv",
        "igphyml_ancestors.common.fa",
        "igphyml_aligned.fa",
        "igphyml_tree.tree",
        "igphyml_tree.pdf"]
    inputs = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/{thing}",
        antibody_lineage = antibody_lineage, chain_type = chain_types, thing = things)
    rule:
        name: f"report_sonar_for_{antibody_lineage}"
        input: inputs
    for chain_type in chain_types:
        inputs = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/{thing}",
            antibody_lineage = antibody_lineage, chain_type = chain_type, thing = things)
        rule:
            name: f"report_sonar_for_{antibody_lineage}.{chain_type}"
            input: inputs

def sonar_island_summary(fp_output_csv, fps_input_csv):
    fieldnames = [
        "specimen", "timepoint", "total", "has_n",
        "germ_div_min", "germ_div_max", "germ_div_median",
        "duplicate_count_median", "cluster_count_median",
        "ab_id_min", "ab_id_max", "ab_id_median"]
    with open(fp_output_csv, "wt") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        rows_out = []
        for fp_in in fps_input_csv:
            rows_out.append(_sonar_island_summary_row(fp_in))
        rows_out = [row for row in rows_out if row["total"] > 0]
        def sorter(row):
            week = re.search("WK([0-9]+)", row["specimen"])
            if week:
                week = int(week.group(1))
            else:
                week = -1
            return (week, row["specimen"])
        rows_out = sorted(rows_out, key=sorter)
        writer.writerows(rows_out)

def _sonar_island_summary_row(fp_in):
    germ_divs = []
    ab_ids_meds = []
    ab_ids_mins = []
    ab_ids_maxes = []
    cluster_counts = []
    duplicate_counts = []
    has_n = 0
    specimen = ""
    timepoint = ""
    with open(fp_in) as f_in:
        reader = csv.DictReader(f_in)
        for row in reader:
            germ_divs.append(float(row["germ_div"]))
            if int(row["n_count"]) > 0:
                has_n += 1
            ab_ids_meds.append(float(row["ab_id_median"]))
            ab_ids_mins.append(float(row["ab_id_min"]))
            ab_ids_maxes.append(float(row["ab_id_max"]))
            duplicate_counts.append(int(row["duplicate_count"]))
            cluster_counts.append(int(row["cluster_count"]))
            # just take the last specimen and timepoint given (if any) since
            # they should be constant per file
            specimen = row.get("specimen", "")
            timepoint = row.get("timepoint", "")
    if germ_divs:
        row_out = {
            "specimen": specimen,
            "timepoint": timepoint,
            "total": len(germ_divs),
            "has_n": has_n,
            "germ_div_min": round(min(germ_divs), 4),
            "germ_div_max": round(max(germ_divs), 4),
            "germ_div_median": round(median(germ_divs), 4),
            "duplicate_count_median": round(median(duplicate_counts)),
            "cluster_count_median": round(median(cluster_counts)),
            "ab_id_min": round(min(ab_ids_mins), 4),
            "ab_id_max": round(max(ab_ids_maxes), 4),
            "ab_id_median": round(median(ab_ids_meds), 4)}
    else:
        # Missing keys will get blanks in the output so that will take care of
        # the rest
        row_out = {
            "specimen": specimen,
            "timepoint": timepoint,
            "total": 0,
            "has_n": 0}
    return row_out

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
    if "specimen" in parts:
        specimens = parts["specimen"]
    else:
        for samp in SAMPLES.values():
            if samp["Type"] == w.chain_type and \
                "IgG" in samp["SpecimenAttrs"]["CellType"]:
                if samp["SpecimenAttrs"]["Subject"] == subject:
                    specimens.add(samp["Specimen"])

    parts["subject"] = subject
    parts["specimen"] = specimens
    # I swear vars(w) *used* to just give you a dictionary of wildcard names
    # and values, but now I'm getting a bunch of functions (and other stuff
    # like _names) mixed in too, which crashes expand().  This is hacky but
    # fixes this for now.
    parts = {key: parts[key] for key in parts if not callable(parts[key])}
    return expand(pattern, **parts)


rule report_sonar_island_summary:
    """Further condense SONAR ID/DIV stats to one file per lineage."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/island_stats_summary.csv"
    input: lambda w: input_helper_sonar(w, "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/{specimen}.island_stats.csv")
    run: sonar_island_summary(output[0], input)

rule report_sonar_island_stats:
    """Condense the full SONAR ID/DIV stats to just those for one island and sumamrize across antibodies."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/{specimen}.island_stats.csv"
    input:
        iddiv=lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.alt.tab"),
        fasta=lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/sequences/nucleotide/islandSeqs_{antibody_lineage}.fa")
    run:
        # just calculate relative to the members of this lineage (for cases
        # where there's more than one)
        mabs = [attrs["Isolate"] for attrs in ANTIBODY_ISOLATES.values() if attrs["Lineage"] == wildcards.antibody_lineage]
        fp_output = output[0]
        fp_input_iddiv = input.iddiv[0]
        fp_input_fasta = input.fasta[0]
        timepoint = [attrs["Timepoint"] for spec, attrs in SPECIMENS.items() if spec == wildcards.specimen][0]
        fieldnames = ["specimen", "timepoint", "sequence_id", "length", "n_count", "v_gene", "germ_div", "ab_id_min", "ab_id_median", "ab_id_max"]
        # outer dict: seq ID to attributes
        # each inner dict: key/val pairs from sequence descriptions
        descs = {}
        seqs = {}
        with open(fp_input_fasta) as f_in:
            for record in SeqIO.parse(f_in, "fasta"):
                descs[record.id] = parse_seq_desc(record.description)
                seqs[record.id] = str(record.seq)
        desc_keys = [val.keys() for val in descs.values()]
        desc_keys = sorted(list(set(itertools.chain(*desc_keys))))
        fieldnames += desc_keys
        with open(fp_input_iddiv) as f_in, open(fp_output, "wt") as f_out:
            reader = csv.DictReader(f_in, delimiter="\t")
            writer = csv.DictWriter(f_out, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            for row in reader:
                if row["sequence_id"] not in descs.keys():
                    continue
                vals = [float(val) for key, val in row.items() if key in mabs]
                if vals:
                    ab_min = min(vals)
                    ab_med = round(median(vals), 4)
                    ab_max = max(vals)
                else:
                    ab_min = ''
                    ab_med = ''
                    ab_max = ''
                row_out = {
                    "specimen": wildcards.specimen,
                    "timepoint": timepoint,
                    "sequence_id": row["sequence_id"],
                    "length": len(seqs[row["sequence_id"]]),
                    "n_count": len(re.sub("[^N]", "", seqs[row["sequence_id"]])),
                    "v_gene": row["v_gene"],
                    "germ_div": row["germ_div"],
                    "ab_id_min": ab_min,
                    "ab_id_median": ab_med,
                    "ab_id_max": ab_max}
                for key in desc_keys:
                    row_out[key] = descs[row["sequence_id"]].get(key, "")
                writer.writerow(row_out)

rule report_sonar_member_rearrangements_summary:
    """Further condense SONAR rearrangements tables to one file per lineage."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/rearrangements_summary.tsv"
    input: lambda w: input_helper_sonar(w, "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/{specimen}.rearrangements.tsv")
    run:
        rows_out = []
        for path in input:
            with open(path) as f_in:
                for row in DictReader(f_in, delimiter="\t"):
                    rows_out.append(row)
        rows_out = sorted(
            rows_out, key=lambda r: (r["timepoint"], r["specimen"], r["sequence_id"]))
        with open(output[0], "wt") as f_out:
            writer = DictWriter(
                f_out, fieldnames=rows_out[0].keys(),
                delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows_out)

rule report_sonar_member_rearrangements:
    """Filter SONAR rearrangements table per specimen to just those for lineage members."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/{specimen}.rearrangements.tsv"
    input:
        seqids=lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/islandSeqs_{antibody_lineage}.txt"),
        tsv=lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_rearrangements.tsv")
    run:
        tp_map = {attrs["Specimen"]: attrs["Timepoint"] for attrs in SPECIMENS.values()}
        with open(input.seqids[0]) as f_in:
            seq_ids = [line.strip() for line in f_in]
        with open(input.tsv[0]) as f_in, open(output[0], "wt") as f_out:
            reader = DictReader(f_in, delimiter="\t")
            fields = ["specimen", "timepoint"] + reader.fieldnames
            writer = DictWriter(
                f_out, fieldnames=fields,
                delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                if row["sequence_id"] in seq_ids:
                    row["specimen"] = wildcards.specimen
                    row["timepoint"] = tp_map[wildcards.specimen]
                    writer.writerow(row)

rule report_sonar_members_table:
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/members.csv"
    input:
        lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/sequences/nucleotide/islandSeqs_{antibody_lineage}.fa")
    run:
        rows = []
        for fasta in input:
            parts = Path(fasta).parts
            specimen = SPECIMENS[parts[3]]
            chain_type = wildcards.chain_type
            chain = "light" if chain_type in ["kappa", "lambda"] else "heavy"
            locus = {"gamma": "H", "lambda": "L", "kappa": "K", "mu": "H"}[chain_type]
            with open(fasta) as f_in:
                for record in SeqIO.parse(f_in, "fasta"):
                    rows.append({
                        "LineageMember": "",
                        "Lineage": wildcards.antibody_lineage,
                        "OriginalID": record.id,
                        "FirstOccurrence": "",
                        "Chain": chain,
                        "Timepoint": int(specimen["Timepoint"]),
                        "TimepointLabel": "",
                        "Specimen": specimen["Specimen"],
                        "Member": "T",
                        "Sequence": str(record.seq)})

        timepoints, labels, rows = format_timepoints(rows)
        for label, row in zip(labels, rows):
            orig_id = f"{label}-{row['OriginalID']}"
            member_name = f"{wildcards.antibody_lineage}-{locus}-{orig_id}"
            row["LineageMember"] = member_name
            row["OriginalID"] = orig_id
            row["TimepointLabel"] = label

        seqmap = defaultdict(list)
        for row in rows:
            seqmap[row["Sequence"]].append((row["TimepointLabel"], row["OriginalID"]))
        for row in rows:
            matches = seqmap[row["Sequence"]]
            if len(matches) > 1:
                matches = sorted(matches)
                if matches[0][1] != row["OriginalID"]:
                    row["FirstOccurrence"] = matches[0][1]
        rows = sorted(rows, key=lambda row: (row["Lineage"], row["Chain"], row["Timepoint"], row["TimepointLabel"], row["LineageMember"]))
        with open(output[0], "wt") as f_out:
            writer = DictWriter(f_out, fieldnames=rows[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

rule report_sonar_igphyml_collected_table:
    """Convert SONAR module 3 collected FASTA into a table."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_collected.csv"
    input: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_collected.fa"
    run:
        rows = []
        with open(input[0]) as f_in:
            for record in SeqIO.parse(f_in, "fasta"):
                fields = parse_seq_desc(record.description)
                tp_label = re.sub("-.*", "", record.id)
                tp_match = re.match(r"wk(N?)([0-9]+)\.?[0-9]*", tp_label)
                tp = ""
                if tp_match:
                    tp = int(tp_match.group(2))
                    if tp_match.group(1) == "N":
                        tp = - tp
                attrs = {
                    "timepoint": tp,
                    "timepoint_label": tp_label,
                    "sequence_id": record.id,
                    "sequence": str(record.seq),
                    }
                attrs.update(fields)
                rows.append(attrs)
        keys = [row.keys() for row in rows]
        keys = list(set(itertools.chain(*keys)))
        field_defaults = [
            "timepoint", "timepoint_label", "num_observations", "num_timepoints", "total_observations",
            "persist", "last_timepoint", "sequence_id", "sequence"]
        def fieldsort(k):
            try:
                return field_defaults.index(k)
            except ValueError:
                return len(keys)
        keys = sorted(keys, key=fieldsort)
        with open(output[0], "w") as f_out:
            writer = DictWriter(f_out, fieldnames=keys, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

rule report_sonar_igphyml_collected:
    """Copy SONAR module 3 collected FASTA.

    This information is subtly different from the per-specimen members.csv
    files, because repeated observations of the same sequence between
    timepoints are collapsed down to the earliest observation and the details
    noted in the sequence descriptions.
    """
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_collected.fa"
    input:
        lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/output/sequences/nucleotide/longitudinal-{antibody_lineage}-collected.fa")
    # Using igseq convert here and elsewhere instead of just cp to unwrap any
    # of these that are wrapped
    shell: "igseq convert {input} {output}"

rule report_sonar_igphyml_ancestors_table:
    """Convert inferred ancestor FASTA into a table with detected clade details.

    This can be used as input to the Inferred sheet.
    """
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.csv"
    input: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.fa"
    run:
        locus = {"gamma": "H", "lambda": "L", "kappa": "K", "mu": "H"}[wildcards.chain_type]
        chain = "heavy" if locus == "H" else "light"
        # names of the antibody isolates for this lineage, as ordered in the
        # metadata.   The IgPhyML output uses all caps.
        isolates = {k.upper(): v for k, v in ANTIBODY_ISOLATES.items() if v["Lineage"] == wildcards.antibody_lineage}
        rows = []
        with open(input[0]) as f_in:
            for record in SeqIO.parse(f_in, "fasta"):
                if re.match(r"1;IG[^;]+;", record.id):
                    continue
                # For each inferred ancestor, we'll check which isolates are
                # within the associated clade and generate an ID based on clade
                # membership.
                tree_index, clade_items, tree_suffix = record.id.split(";")
                clade_items = clade_items.split(",")
                mask = sum([(iso not in clade_items)*2**(exp) for exp, iso in enumerate(isolates.keys())])
                # mask uses a 1 bit to mean NOT in the clade, so all 1 means no
                # isolates are present.
                maxmask = 2**len(isolates) - 1
                if mask == maxmask:
                    mask = ""
                else:
                    mask = hex(mask)[2:].upper()
                # Timepoint
                timepoints_isolate = [isolates.get(item, {}).get("Timepoint") for item in clade_items]
                timepoint_min = None
                for item in clade_items:
                    match = re.match(r"WK([0-9]+)(\.?[0-9]*)-[0-9]+", item)
                    if match:
                        timepoint = int(match.group(1))
                        if timepoint_min is None or timepoint < timepoint_min:
                            timepoint_min = timepoint
                    else:
                        timepoint = isolates.get(item, {}).get("Timepoint")
                        if timepoint:
                            timepoint = int(timepoint)
                            if timepoint_min is None or timepoint < timepoint_min:
                                timepoint_min = timepoint
                # Name, for ancestors that lead to any of the isolates
                name = ""
                if mask != "":
                    name_fields = [wildcards.antibody_lineage, locus, tree_index]
                    if mask != "0":
                        name_fields.append(mask)
                    name = "-".join(name_fields)
                # Append to one big list so we can sort it before writing
                rows.append({
                    "InferredAncestor": name,
                    "Lineage": wildcards.antibody_lineage,
                    "TreeDepth": int(tree_index),
                    "Chain": chain,
                    "IsolateSubsetID": mask,
                    "Timepoint": int(timepoint_min),
                    "OriginalID": record.id,
                    "Sequence": str(record.seq)})
        rows = sorted(rows, key=lambda r: (r["Timepoint"], r["TreeDepth"], r["IsolateSubsetID"]))
        with open(output[0], "wt") as f_out:
            fields = [
                "InferredAncestor", "Lineage", "TreeDepth", "Chain",
                "IsolateSubsetID", "Timepoint", "OriginalID", "Sequence"]
            writer = DictWriter(f_out, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

rule report_sonar_igphyml_ancestors_common:
    """Make version of SONAR module 3 inferred ancestors FASTA filtered to mAb ancestors."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.common.fa"
    input:
        ancs="analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.fa",
        mabs="analysis/reporting/sonar/{antibody_lineage}.{chain_type}/mabs.csv"
    shell: "sonar_ancs_common.py --ancestors {input.ancs} --clade {input.mabs} --output {output}"

rule report_sonar_igphyml_ancestors_custom_common:
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.custom.common.fa"
    input:
        ancs="analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.custom.fa",
        mabs="analysis/reporting/sonar/{antibody_lineage}.{chain_type}/mabs.csv"
    shell: "sonar_ancs_common.py --ancestors {input.ancs} --clade {input.mabs} --output {output}"

rule report_sonar_igphyml_ancestors:
    """Copy SONAR module 3 inferred ancestors FASTA."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.fa"
    input:
        lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/output/sequences/nucleotide/longitudinal-{antibody_lineage}_inferredAncestors.fa")
    shell: "igseq convert {input} {output}"

rule report_sonar_igphyml_ancestors_custom:
    """Copy SONAR module 3 (custom alignment) inferred ancestors FASTA."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.custom.fa"
    input:
        lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/longitudinal-custom-{antibody_lineage}/output/sequences/nucleotide/longitudinal-custom-{antibody_lineage}_inferredAncestors.fa")
    shell: "igseq convert {input} {output}"

rule report_sonar_igphyml_alignment:
    """Copy SONAR module 3 alignment FASTA."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_aligned.fa"
    input:
        lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/work/phylo/longitudinal-{antibody_lineage}_aligned.afa")
    shell: "igseq convert {input} {output}"

rule report_sonar_igphyml_tree_pdf:
    """Copy SONAR module 3 tree PDF."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_tree.pdf"
    input:
        lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/output/longitudinal-{antibody_lineage}_igphyml.tree.pdf")
    shell: "cp {input} {output}"

rule report_sonar_igphyml_tree:
    """Copy SONAR module 3 newick tree file."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_tree.tree"
    input:
        lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/output/longitudinal-{antibody_lineage}_igphyml.tree")
    shell: "cp {input} {output}"

rule report_sonar_mabs:
    # not a SONAR output, exactly, but other outputs are relative to this
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/mabs.csv"
    run:
        chain = "heavy"
        if wildcards.chain_type in ["kappa", "lambda"]:
            chain = "light"
        with open(output[0], "w") as f_out:
            writer = DictWriter(
                f_out,
                fieldnames=["timepoint", "chain", "sequence_id", "sequence"],
                lineterminator="\n")
            writer.writeheader()
            for attrs in ANTIBODY_ISOLATES.values():
                if attrs["Lineage"] == wildcards.antibody_lineage:
                    seq = attrs[chain.capitalize() + "Seq"]
                    if seq:
                        writer.writerow({
                            "timepoint": attrs["Timepoint"],
                            "chain": chain,
                            "sequence_id": attrs["Isolate"],
                            "sequence": seq})
