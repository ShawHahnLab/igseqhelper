import itertools
from statistics import median
from collections import defaultdict
from csv import DictWriter, DictReader
from Bio import SeqIO

rule all_counts:
    input: expand("analysis/reporting/counts/counts_by_{thing}.csv", thing=["sample", "run", "specimen"])

# this rule has no inputs so it'll just take whatever files are available at
# run time.
rule counts_by_sample:
    output: "analysis/reporting/counts/counts_by_sample.csv"
    run: counts_by_sample(output[0])

rule counts_by_run:
    output: "analysis/reporting/counts/counts_by_run.csv"
    input: "analysis/reporting/counts/counts_by_sample.csv"
    run: counts_by_run(input[0], output[0])

rule counts_by_specimen:
    output: "analysis/reporting/counts/counts_by_specimen.csv"
    input: "analysis/reporting/counts/counts_by_sample.csv"
    run: counts_by_specimen(input[0], output[0])

def counts_for(samp, runid, category):
    path = Path("analysis")/category/runid/f"{samp}.{category}.counts.csv"
    if path.exists():
        with open(path) as f_in:
            reader = DictReader(f_in)
            for row in reader:
                if row["Item"] == "output":
                    return row["NumSeqs"]
    return ""

def divide(val1, val2, fmt="{:.2f}"):
    """Divide val1 by val2 as floats and return formatted string."""
    try:
        num = float(val1)/float(val2)
    except (ValueError, TypeError, ZeroDivisionError):
        num = ""
    else:
        num = fmt.format(num)
    return num

def counts_by_sample(csv_out):
    fieldnames = [
        "Run", "Specimen", "CellType", "Type", "Sample",
        "CountsDemux", "CountsTrim", "CountsMerge",
        "CellCount", "RatioDemux", "RatioMerge"]
    demux = {}
    phix = {}
    rows_out = []
    with open(csv_out, "wt") as f_out:
        writer = DictWriter(
            f_out,
            fieldnames=fieldnames,
            lineterminator="\n")
        writer.writeheader()
        for samp, attrs in SAMPLES.items():
            path_demux = Path("analysis/demux")/attrs["Run"]/"demux.counts.csv"
            path_phix = Path("analysis/phix")/attrs["Run"]/"phix.counts.csv"
            samp_demux = ""
            if path_demux.exists():
                if attrs["Run"] not in demux:
                    with open(path_demux) as f_in:
                        rows = list(DictReader(f_in))
                    demux[attrs["Run"]] = rows
                for row in demux[attrs["Run"]]:
                    if row.get("Sample") == samp:
                        samp_demux = row["NumSeqs"]
            if path_phix.exists():
                if attrs["Run"] not in phix:
                    with open(path_phix) as f_in:
                        rows = list(DictReader(f_in))
                    phix[attrs["Run"]] = rows
            samp_trim = counts_for(samp, attrs["Run"], "trim")
            samp_merge = counts_for(samp, attrs["Run"], "merge")
            row = {
                "Run": attrs["Run"],
                "Specimen": attrs["Specimen"],
                "CellType": attrs["SpecimenAttrs"]["CellType"],
                "Type": attrs["Type"],
                "Sample": samp,
                "CountsDemux": samp_demux,
                "CountsTrim": samp_trim,
                "CountsMerge": samp_merge,
                "CellCount": attrs["SpecimenAttrs"]["CellCount"],
                "RatioDemux": divide(samp_demux, attrs["SpecimenAttrs"]["CellCount"]),
                "RatioMerge": divide(samp_merge, attrs["SpecimenAttrs"]["CellCount"])
                }
            rows_out.append(row)
        # write per-run info
        for runid, rows in demux.items():
            for row in rows:
                if row.get("Item") == "unassigned":
                    rows_out.append({"Run": runid, "Sample": "unassigned", "CountsDemux": row["NumSeqs"]})
        for runid, rows in phix.items():
            for row in rows:
                if row.get("Item") == "mapped":
                    rows_out.append({"Run": runid, "Sample": "unassigned.phix", "CountsDemux": row["NumSeqs"]})
        # TODO phix
        rows_out = sorted(rows_out, key=lambda r: (str(r.get("Run")), str(r.get("Sample"))))
        writer.writerows(rows_out)

def counts_by_run(input_csv, output_csv):
    run_info = {}
    with open(input_csv) as f_in:
        reader = DictReader(f_in)
        for row in reader:
            if row["Run"]:
                if row["Run"] not in run_info:
                    run_info[row["Run"]] = {
                        "Run": row["Run"],
                        "UnassignedSeqs": "",
                        "PhixSeqs": "",
                        "SampleSeqs": []}
                path_run_counts = Path("analysis/reads")/row["Run"]/"getreads.counts.csv"
                if path_run_counts.exists():
                    with open(path_run_counts) as f_in:
                        for countrow in DictReader(f_in):
                            if countrow["Item"] == "unassigned-raw":
                                run_info[row["Run"]]["RawReads"] = countrow["NumSeqs"]
                            elif countrow["Item"] == "unassigned-pf":
                                run_info[row["Run"]]["PassingFilter"] = countrow["NumSeqs"]
                if row["Sample"] == "unassigned":
                    if row["CountsDemux"]:
                        run_info[row["Run"]]["UnassignedSeqs"] = int(row["CountsDemux"])
                elif row["Sample"] == "unassigned.phix":
                    if row["CountsDemux"]:
                        run_info[row["Run"]]["PhixSeqs"] = int(row["CountsDemux"])
                else:
                    if row["CountsDemux"]:
                        run_info[row["Run"]]["SampleSeqs"].append(int(row["CountsDemux"]))
    with open(output_csv, "wt") as f_out:
        writer = DictWriter(
            f_out,
            fieldnames=["Run", "RawReads", "PassingFilter", "UnassignedSeqs", "PhixSeqs", "SampleSeqs", "UnassignedFraction", "PhixFraction"],
            lineterminator="\n")
        writer.writeheader()
        for row in run_info.values():
            row["SampleSeqs"] = sum(row["SampleSeqs"]) if row["SampleSeqs"] else ""
            parts = [row["UnassignedSeqs"], row["SampleSeqs"]]
            parts = [p for p in parts if p]
            if parts:
                row["TotalSeqs"] = sum(parts)
            row["UnassignedFraction"] = divide(row["UnassignedSeqs"], row["TotalSeqs"])
            row["PhixFraction"] = divide(row["PhixSeqs"], row["TotalSeqs"])
            del row["TotalSeqs"]
            writer.writerow(row)

def counts_by_specimen(input_csv, output_csv):
    spec_info = {}
    with open(input_csv) as f_in:
        reader = DictReader(f_in)
        for row in reader:
            key = (row["Specimen"], row["Type"])
            if row["Specimen"]:
                if key not in spec_info:
                    spec_info[key] = {
                        "Specimen": row["Specimen"],
                        "CellType": row["CellType"],
                        "CellCount": row["CellCount"],
                        "Type": row["Type"],
                        "DemuxSeqs": [], "TrimSeqs": [], "MergeSeqs": []}
                if row["CountsDemux"]:
                    spec_info[key]["DemuxSeqs"].append(int(row["CountsDemux"]))
                if row["CountsTrim"]:
                    spec_info[key]["TrimSeqs"].append(int(row["CountsTrim"]))
                if row["CountsMerge"]:
                    spec_info[key]["MergeSeqs"].append(int(row["CountsMerge"]))
    with open(output_csv, "wt") as f_out:
        writer = DictWriter(
            f_out,
            fieldnames=["Specimen", "CellType", "Type", "DemuxSeqs", "TrimSeqs", "MergeSeqs", "CellCount", "RatioDemux", "RatioMerge"],
            lineterminator="\n")
        writer.writeheader()
        for row in spec_info.values():
            row["DemuxSeqs"] = sum(row["DemuxSeqs"]) if row["DemuxSeqs"] else ""
            row["TrimSeqs"] = sum(row["TrimSeqs"]) if row["TrimSeqs"] else ""
            row["MergeSeqs"] = sum(row["MergeSeqs"]) if row["MergeSeqs"] else ""
            row["RatioDemux"] = divide(row["DemuxSeqs"], row["CellCount"])
            row["RatioMerge"] = divide(row["MergeSeqs"], row["CellCount"])
            writer.writerow(row)

### SONAR Members and Ancestors

rule sonar_island_stats:
    """Condense the full ID/DIV stats to just those for one island and sumamrize across antibodies."""
    output: "analysis/reporting/by-lineage/{antibody_lineage}.{chain_type}/{specimen}.island_stats.csv"
    input:
        iddiv=lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab"),
        fasta=lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/{specimen}/output/sequences/nucleotide/islandSeqs_{antibody_lineage}.fa")
    run:
        # just calculate relative to the members of this lineage (for cases
        # where there's more than one)
        mabs = [attrs["AntibodyIsolate"] for attrs in ANTIBODY_ISOLATES.values() if attrs["AntibodyLineage"] == wildcards.antibody_lineage]
        fp_output = output[0]
        fp_input_iddiv = input.iddiv[0]
        fp_input_fasta = input.fasta[0]
        timepoint = [attrs["Timepoint"] for spec, attrs in SPECIMENS.items() if spec == wildcards.specimen][0]
        fieldnames = ["specimen", "timepoint", "sequence_id", "v_gene", "germ_div", "ab_id_min", "ab_id_median", "ab_id_max"]
        # outer dict: seq ID to attributes
        # each inner dict: key/val pairs from sequence descriptions
        with open(fp_input_fasta) as f_in:
            descs = {rec.id: igseqhelper.util.parse_seq_desc(rec.description) for rec in SeqIO.parse(f_in, "fasta")}
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
                    "v_gene": row["v_gene"],
                    "germ_div": row["germ_div"],
                    "ab_id_min": ab_min,
                    "ab_id_median": ab_med,
                    "ab_id_max": ab_max}
                for key in desc_keys:
                    row_out[key] = descs[row["sequence_id"]].get(key, "")
                writer.writerow(row_out)

def sonar_island_summary(fp_output_csv, fps_input_csv):
    fieldnames = [
        "specimen", "timepoint", "total",
        "germ_div_min", "germ_div_max", "germ_div_median",
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
    specimen = ""
    timepoint = ""
    with open(fp_in) as f_in:
        reader = csv.DictReader(f_in)
        for row in reader:
            germ_divs.append(float(row["germ_div"]))
            ab_ids_meds.append(float(row["ab_id_median"]))
            ab_ids_mins.append(float(row["ab_id_min"]))
            ab_ids_maxes.append(float(row["ab_id_max"]))
            # just take the last specimen and timepoint given (if any) since
            # they should be constant per file
            specimen = row.get("specimen", "")
            timepoint = row.get("timepoint", "")
    if germ_divs:
        row_out = {
            "specimen": specimen,
            "timepoint": timepoint,
            "total": len(germ_divs),
            "germ_div_min": round(min(germ_divs), 4),
            "germ_div_max": round(max(germ_divs), 4),
            "germ_div_median": round(median(germ_divs), 4),
            "ab_id_min": round(min(ab_ids_mins), 4),
            "ab_id_max": round(max(ab_ids_maxes), 4),
            "ab_id_median": round(median(ab_ids_meds), 4)}
    else:
        row_out = {
            "specimen": specimen,
            "timepoint": timepoint,
            "total": 0,
            "germ_div_min": '',
            "germ_div_max": '',
            "germ_div_median": '',
            "ab_id_min": '',
            "ab_id_max": '',
            "ab_id_median": ''}
    return row_out

rule sonar_island_summary:
    """Further condense ID/DIV stats to one file per lineage."""
    output: "analysis/reporting/by-lineage/{antibody_lineage}.{chain_type}/island_stats_summary.csv"
    input: lambda w: input_helper_sonar(w, "analysis/reporting/by-lineage/{antibody_lineage}.{chain_type}/{specimen}.island_stats.csv")
    run: sonar_island_summary(output[0], input)

rule sonar_members_table:
    output: "analysis/reporting/by-lineage/{antibody_lineage}.{chain_type}/members.csv"
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
            timepoint = specimen["Timepoint"]
            with open(fasta) as f_in:
                for record in SeqIO.parse(f_in, "fasta"):
                    orig_id = f"wk{timepoint}-{record.id}"
                    rows.append({
                        "LineageMember": f"{wildcards.antibody_lineage}-{locus}-{orig_id}",
                        "AntibodyLineage": wildcards.antibody_lineage,
                        "OriginalID": orig_id,
                        "FirstOccurrence": "",
                        "Chain": chain,
                        "Timepoint": int(timepoint),
                        "Member": "T",
                        "Sequence": str(record.seq)})
        seqmap = defaultdict(list)
        for row in rows:
            seqmap[row["Sequence"]].append((row["Timepoint"], row["LineageMember"]))
        for row in rows:
            matches = seqmap[row["Sequence"]]
            if len(matches) > 1:
                matches = sorted(matches)
                row["FirstOccurrence"] = matches[0][1]
        rows = sorted(rows, key=lambda row: (row["AntibodyLineage"], row["Chain"], row["Timepoint"], row["LineageMember"]))
        with open(output[0], "wt") as f_out:
            writer = DictWriter(f_out, fieldnames=rows[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


rule sonar_ancestors_table:
    """Convert inferred ancestor FASTA into a table with detected clade details.

    This can be used as input to the Inferred sheet.
    """
    output: "analysis/reporting/by-lineage/{antibody_lineage}.{chain_type}/ancestors.csv"
    input: lambda w: expand(\
        "analysis/sonar/{subject}.{{chain_type}}/longitudinal-{{antibody_lineage}}/output/sequences/nucleotide/longitudinal-{{antibody_lineage}}_inferredAncestors.fa", \
        subject=ANTIBODY_LINEAGES[w.antibody_lineage]["Subject"])
    run:
        locus = {"gamma": "H", "lambda": "L", "kappa": "K", "mu": "H"}[wildcards.chain_type]
        chain = "heavy" if locus == "H" else "light"
        # names of the antibody isolates for this lineage, as ordered in the
        # metadata.   The IgPhyML output uses all caps.
        isolates = {k.upper(): v for k, v in ANTIBODY_ISOLATES.items() if v["AntibodyLineage"] == wildcards.antibody_lineage}
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
                    "AntibodyLineage": wildcards.antibody_lineage,
                    "TreeDepth": int(tree_index),
                    "Chain": chain,
                    "IsolateSubsetID": mask,
                    "Timepoint": int(timepoint_min),
                    "OriginalID": record.id,
                    "Sequence": str(record.seq)})
        rows = sorted(rows, key=lambda r: (r["Timepoint"], r["TreeDepth"], r["IsolateSubsetID"]))
        with open(output[0], "wt") as f_out:
            fields = [
                "InferredAncestor", "AntibodyLineage", "TreeDepth", "Chain",
                "IsolateSubsetID", "Timepoint", "OriginalID", "Sequence"]
            writer = DictWriter(f_out, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
