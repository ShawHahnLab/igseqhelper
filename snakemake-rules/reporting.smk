import itertools
from statistics import median
from collections import defaultdict
from csv import DictWriter, DictReader
from Bio import SeqIO
import igseqhelper.util

# Just a helper for the below rules.
# Filter out any entries that don't match the metadata since we rely on the
# metadata when prepping the tables.  (This could be things like
# "subject-(datestamp)" or what have you)
def filter_wildcards_by_metadata(wildcards_dict):
    removes = set()
    for key, vec in wildcards_dict.items():
        for idx, item in enumerate(vec):
            if (key == "subject" and item not in {row["Subject"] for row in SPECIMENS.values()}) or \
                (key == "chain_type" and item not in {row["Type"] for row in SAMPLES.values()}) or \
                (key == "specimen" and item not in SPECIMENS.keys()) or \
                (key == "antibody_lineage" and item not in ANTIBODY_LINEAGES.keys()):
                removes.add(idx)
    # (careful, need to remove consistently across all wildcards)
    wildcards_dict_out = {}
    for key, vec in wildcards_dict.items():
        wildcards_dict_out[key] = [item for idx, item in enumerate(vec) if idx not in removes]
    return wildcards_dict_out

AVAILABLE_SONAR_ISLANDS = filter_wildcards_by_metadata(dict(zip(
    ["subject", "chain_type", "specimen", "antibody_lineage", "dummy"],
    glob_wildcards("analysis/sonar/{subject}.{chain_type}/{specimen}/"
        "output/sequences/nucleotide/islandSeqs_{dummy}.fa"))))

AVAILABLE_SONAR_ANCESTORS = filter_wildcards_by_metadata(dict(zip(
    ["subject", "chain_type", "antibody_lineage", "dummy"],
    glob_wildcards("analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/"
        "output/sequences/nucleotide/longitudinal-{dummy}_inferredAncestors.fa"))))

AVAILABLE_IGDISCOVER = dict(zip(
    ["ref", "chain_type", "subject", "segment"],
    glob_wildcards("analysis/igdiscover/{ref}/{chain_type}/{subject}/final/database/{segment}.fasta")))

AVAILABLE_MININGD = filter_wildcards_by_metadata(dict(zip(
    ["subject"], glob_wildcards("analysis/mining-d/{subject}.output.fasta"))))

TARGET_REPORT_COUNTS = expand("analysis/reporting/counts/counts_by_{thing}.csv", thing=["sample", "run", "specimen"])
TARGET_REPORT_SONAR_ISLAND_SUMMARIES = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/island_stats_summary.csv", zip, **AVAILABLE_SONAR_ISLANDS)
TARGET_REPORT_SONAR_MEMBERS_TABLES = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/members.csv", zip, **AVAILABLE_SONAR_ISLANDS)
TARGET_REPORT_SONAR_IGPHYML_COLLECTED = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_collected.csv", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_SONAR_IGPHYML_ANCESTORS = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.csv", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_SONAR_IGPHYML_ANCS_COMMON = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.common.csv", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_SONAR_IGPHYML_ALIGNED = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_aligned.fa", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_SONAR_IGPHYML_TREES = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_tree.tree", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_SONAR_IGPHYML_TREEPDFS = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_tree.pdf", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_IGDISCOVER_TREES = expand("analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.nex", zip, **AVAILABLE_IGDISCOVER)
TARGET_REPORT_IGDISCOVER_FINAL_DBS = expand("analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.fasta", zip, **AVAILABLE_IGDISCOVER)
TARGET_REPORT_MININGD_FASTAS = expand("analysis/reporting/mining-d/{subject}/{subject}.fasta", zip, **AVAILABLE_MININGD)
TARGET_REPORT_MININGD_TREES = expand("analysis/reporting/mining-d/{subject}/{subject}.nex", zip, **AVAILABLE_MININGD)

rule all_report:
    input: TARGET_REPORT_COUNTS +
        TARGET_REPORT_SONAR_ISLAND_SUMMARIES +
        TARGET_REPORT_SONAR_MEMBERS_TABLES +
        TARGET_REPORT_SONAR_IGPHYML_COLLECTED +
        TARGET_REPORT_SONAR_IGPHYML_ANCESTORS +
        TARGET_REPORT_SONAR_IGPHYML_ANCS_COMMON +
        TARGET_REPORT_SONAR_IGPHYML_ALIGNED +
        TARGET_REPORT_SONAR_IGPHYML_TREES +
        TARGET_REPORT_SONAR_IGPHYML_TREEPDFS +
        TARGET_REPORT_IGDISCOVER_TREES +
        TARGET_REPORT_IGDISCOVER_FINAL_DBS +
        TARGET_REPORT_MININGD_TREES +
        ["analysis/reporting/mining-d/all.nex"]

rule report_available_sonar:
    input: TARGET_REPORT_SONAR_ISLAND_SUMMARIES +
        TARGET_REPORT_SONAR_MEMBERS_TABLES +
        TARGET_REPORT_SONAR_IGPHYML_COLLECTED +
        TARGET_REPORT_SONAR_IGPHYML_ANCESTORS +
        TARGET_REPORT_SONAR_IGPHYML_ANCS_COMMON +
        TARGET_REPORT_SONAR_IGPHYML_ALIGNED +
        TARGET_REPORT_SONAR_IGPHYML_TREES +
        TARGET_REPORT_SONAR_IGPHYML_TREEPDFS

rule report_available_igdiscover:
    input: TARGET_REPORT_IGDISCOVER_TREES +
        TARGET_REPORT_IGDISCOVER_FINAL_DBS

rule report_counts:
    input: TARGET_REPORT_COUNTS

rule report_available_sonar_island_summaries:
    input: TARGET_REPORT_SONAR_ISLAND_SUMMARIES

rule report_available_sonar_members_tables:
    input: TARGET_REPORT_SONAR_MEMBERS_TABLES

rule report_available_sonar_igphyml_collected:
    input: TARGET_REPORT_SONAR_IGPHYML_COLLECTED

rule report_available_sonar_igphyml_ancestors:
    input: TARGET_REPORT_SONAR_IGPHYML_ANCESTORS

rule report_available_sonar_igphyml_ancs_common:
    input: TARGET_REPORT_SONAR_IGPHYML_ANCS_COMMON

rule report_available_sonar_igphyml_aligned:
    input: TARGET_REPORT_SONAR_IGPHYML_ALIGNED

rule report_available_sonar_igphyml_trees:
    input: TARGET_REPORT_SONAR_IGPHYML_TREES

rule report_available_sonar_igphyml_treepdfs:
    input: TARGET_REPORT_SONAR_IGPHYML_TREEPDFS

rule report_available_igdiscover_trees:
    input: TARGET_REPORT_IGDISCOVER_TREES

rule report_available_igdiscover_final_dbs:
    input: TARGET_REPORT_IGDISCOVER_FINAL_DBS

rule report_available_miningd_fastas:
    input: TARGET_REPORT_MININGD_FASTAS

rule report_available_miningd_trees:
    input: TARGET_REPORT_MININGD_TREES

### Counts

# this rule has no inputs so it'll just take whatever files are available at
# run time.
rule report_counts_by_sample:
    output: "analysis/reporting/counts/counts_by_sample.csv"
    run: counts_by_sample(output[0])

rule report_counts_by_run:
    output: "analysis/reporting/counts/counts_by_run.csv"
    input: "analysis/reporting/counts/counts_by_sample.csv"
    run: counts_by_run(input[0], output[0])

rule report_counts_by_specimen:
    output: "analysis/reporting/counts/counts_by_specimen.csv"
    input: "analysis/reporting/counts/counts_by_sample.csv"
    run: counts_by_specimen(input[0], output[0])

def report_counts_for(samp, runid, category):
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
        "Run", "Subject", "Specimen", "CellType", "Type", "Sample",
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
            samp_trim = report_counts_for(samp, attrs["Run"], "trim")
            samp_merge = report_counts_for(samp, attrs["Run"], "merge")
            row = {
                "Run": attrs["Run"],
                "Subject": attrs["SpecimenAttrs"]["Subject"],
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
            row["TotalSeqs"] = sum(parts) if parts else ""
            row["UnassignedFraction"] = divide(row["UnassignedSeqs"], row["TotalSeqs"])
            row["PhixFraction"] = divide(row["PhixSeqs"], row["TotalSeqs"])
            del row["TotalSeqs"]
            writer.writerow(row)

def sonar_airr_counts(input_tsv, output_csv, fmt_islands=None, lineages=None):
    # Reads: grand total in the file.  Some don't make it this far if they
    #        didn't even look like antibody reads (e.g. ferritin)
    # GoodReads: The good (no stops, all parts intact) antibody reads
    # ClusteredReads: Good ones that made it into clusters (not quite all do)
    # ClusteredUnique The *number* of clusters
    # LineageMembers: sum of lineage members across any lineages for this
    #                 specimen (only included if files are available)
    counts = {"SONARReads": 0, "SONARGoodReads": 0, "SONARClusteredReads": 0, "SONARClusteredUnique": 0}
    with open(input_tsv) as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        for row in reader:
            counts["SONARReads"] += int(row["duplicate_count"])
            if row["status"] == "good":
                counts["SONARGoodReads"] += int(row["duplicate_count"])
                if row["cluster_count"]:
                    # counting the number of reads counted in clusters.  Ony good
                    # reads get clustered, and only those with the minimum dupliate
                    # count (by default singletons are not assigned clusters, good
                    # or not).
                    counts["SONARClusteredReads"] += int(row["cluster_count"])
                    # counting each cluster once for this one
                    counts["SONARClusteredUnique"] += 1
    # If there are islandSeqs files for the expected lineages, tally those up
    # too (but just if available)
    if lineages and fmt_islands:
        members = None
        for antibody_lineage in lineages:
            path = Path(fmt_islands.format(antibody_lineage = antibody_lineage))
            if path.exists():
                members = members or 0
                with open(path) as f_in:
                    for line in f_in:
                        members += 1
        if members is not None:
            counts["LineageMembers"] = members
    with open(output_csv, "wt") as f_out:
        writer = DictWriter(f_out, fieldnames=counts.keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerow(counts)

rule report_sonar_airr_counts:
    output: csv="analysis/reporting/counts/sonar/{subject}.{chain_type}.{specimen}.csv"
    input: tsv="analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_rearrangements.tsv"
    run:
        fmt_islands = "analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/islandSeqs_{{antibody_lineage}}.txt"
        fmt_islands = fmt_islands.format(**vars(wildcards))
        lineages = [lineage for lineage, attrs in ANTIBODY_LINEAGES.items() if attrs["Subject"] == wildcards.subject]
        sonar_airr_counts(input.tsv, output.csv, fmt_islands, lineages)

rule report_sonar_airr_counts_by_subject_and_type:
    output: touch("analysis/reporting/counts/sonar/{subject}.{chain_type}.done")
    input: lambda w: expand("analysis/reporting/counts/sonar/{{subject}}.{{chain_type}}.{specimen}.csv", specimen={samp["Specimen"] for samp in SAMPLES.values() if samp["SpecimenAttrs"]["Subject"] == w.subject and samp["Type"] == w.chain_type})

def counts_by_specimen(input_csv, output_csv):
    spec_info = {}
    with open(input_csv) as f_in:
        reader = DictReader(f_in)
        for row in reader:
            key = (row["Subject"], row["Specimen"], row["Type"])
            if row["Specimen"]:
                if key not in spec_info:
                    spec_info[key] = {
                        "Subject": row["Subject"],
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
                # These take a while to crunch through so I'm doing that in a
                # separate rule, but still not explicitly giving it as input so
                # this rule will only use whatever's already on disk
                sonar_counts = Path(f"analysis/reporting/counts/sonar/{row['Subject']}.{row['Type']}.{row['Specimen']}.csv")
                if sonar_counts.exists():
                    with open(sonar_counts) as f_in_sonar:
                        reader = DictReader(f_in_sonar)
                        spec_info[key].update(next(reader))
    fieldnames = ["Subject", "Specimen", "CellType", "Type",
        "DemuxSeqs", "TrimSeqs", "MergeSeqs", "CellCount", "RatioDemux", "RatioMerge",
        "SONARReads", "SONARGoodReads", "SONARClusteredReads", "SONARClusteredUnique",
        "LineageMembers"]
    rows = []
    for row in spec_info.values():
        row["DemuxSeqs"] = sum(row["DemuxSeqs"]) if row["DemuxSeqs"] else ""
        row["TrimSeqs"] = sum(row["TrimSeqs"]) if row["TrimSeqs"] else ""
        row["MergeSeqs"] = sum(row["MergeSeqs"]) if row["MergeSeqs"] else ""
        row["RatioDemux"] = divide(row["DemuxSeqs"], row["CellCount"])
        row["RatioMerge"] = divide(row["MergeSeqs"], row["CellCount"])
        rows.append(row)
    rows = sorted(rows, key=lambda row: [row.get(field, "") for field in fieldnames])
    with open(output_csv, "wt") as f_out:
        writer = DictWriter(
            f_out,
            fieldnames=fieldnames,
            lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

### Lineages

def report_lineages_divergence_input(w):
    # queries
    chain_types = ["gamma"]
    for attrs in ANTIBODY_LINEAGES.values():
        if attrs["AntibodyLineage"] == w.antibody_lineage:
            try:
                chain_types.append({"L": "lambda", "K": "kappa"}[attrs["VL"][2]])
            except:
                pass
    targets = {
        "members": expand(
            "analysis/reporting/sonar/{{antibody_lineage}}.{chain_type}/igphyml_collected.csv",
            chain_type = chain_types),
        "mabs": ["analysis/reporting/by-lineage/{antibody_lineage}.mabs.csv"]}
    # refs
    mktargets = lambda ct, segs: expand(
        "analysis/reporting/igdiscover/sonarramesh/{chain_type}/{subject}/{segment}.fasta",
        antibody_lineage=w.antibody_lineage, chain_type=ct, subject=subject, segment=segs)
    light_ct = None
    for attrs in ANTIBODY_LINEAGES.values():
        if attrs["AntibodyLineage"] == w.antibody_lineage:
            subject = attrs["Subject"]
            try:
                light_ct = {"L": "lambda", "K": "kappa"}[attrs["VL"][2]]
                break
            except:
                pass
    else:
        raise ValueError
    targets["refs"] = mktargets("mu", ["V", "D", "J"])
    if light_ct:
        targets["refs"] += mktargets(light_ct, ["V", "J"])
    return targets

def report_lineages_divergence_param_refs(w, input):
    # to condense V/D/J to just parent dirs
    return list({Path(path).parent for path in input.refs})

rule report_lineages_divergence_plot:
    output: "analysis/reporting/by-lineage/{antibody_lineage}.divergence.pdf"
    input: "analysis/reporting/by-lineage/{antibody_lineage}.divergence.csv"
    shell: "germ_div.py -Q {input} -o {output}"

rule report_lineages_divergence:
    output: "analysis/reporting/by-lineage/{antibody_lineage}.divergence.csv"
    input: unpack(report_lineages_divergence_input)
    params:
        refs=report_lineages_divergence_param_refs
    shell: "germ_div.py -S rhesus -r {params.refs} -Q {input.members} mabs={input.mabs} -G Member mabs=mAb -o {output}"

rule report_lineages_mabs:
    output: temp("analysis/reporting/by-lineage/{antibody_lineage}.mabs.csv")
    run:
        with open(output[0], "w") as f_out:
            writer = DictWriter(
                f_out,
                fieldnames=["timepoint", "chain", "sequence_id", "sequence"],
                lineterminator="\n")
            writer.writeheader()
            for attrs in ANTIBODY_ISOLATES.values():
                if attrs["AntibodyLineage"] == wildcards.antibody_lineage:
                    for chain in ["heavy", "light"]:
                        seq = attrs[chain.capitalize() + "Seq"]
                        if seq:
                            writer.writerow({
                                "timepoint": attrs["Timepoint"],
                                "chain": chain,
                                "sequence_id": attrs["AntibodyIsolate"] + f"_{chain}",
                                "sequence": seq})

### SONAR Members and Ancestors

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
        mabs = [attrs["AntibodyIsolate"] for attrs in ANTIBODY_ISOLATES.values() if attrs["AntibodyLineage"] == wildcards.antibody_lineage]
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
                descs[record.id] = igseqhelper.util.parse_seq_desc(record.description)
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
                        "AntibodyLineage": wildcards.antibody_lineage,
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
        rows = sorted(rows, key=lambda row: (row["AntibodyLineage"], row["Chain"], row["Timepoint"], row["TimepointLabel"], row["LineageMember"]))
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
                fields = igseqhelper.util.parse_seq_desc(record.description)
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

rule report_sonar_igphyml_ancestors_common:
    """Make verison of SONAR module 3 inferred ancestors FASTA filtered to mAb ancestors."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.common.fa"
    input:
        ancs="analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.fa",
        mabs="analysis/reporting/sonar/{antibody_lineage}.{chain_type}/mabs.csv"
    shell: "sonar_ancs_common.py --ancestors {input.ancs} --clade {input.mabs} --output {output}"

rule report_sonar_igphyml_ancestors:
    """Copy SONAR module 3 inferred ancestors FASTA."""
    output: "analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.fa"
    input:
        lambda w: input_helper_sonar(w, "analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/output/sequences/nucleotide/longitudinal-{antibody_lineage}_inferredAncestors.fa")
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
                if attrs["AntibodyLineage"] == wildcards.antibody_lineage:
                    seq = attrs[chain.capitalize() + "Seq"]
                    if seq:
                        writer.writerow({
                            "timepoint": attrs["Timepoint"],
                            "chain": chain,
                            "sequence_id": attrs["AntibodyIsolate"],
                            "sequence": seq})

### IgDiscover

# See also the final/dendrogram_{segment}.pdf files from IgDiscover
rule report_igdiscover_tree:
    output: "analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.nex"
    input:
        after="analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.fasta",
        before="analysis/igdiscover/{ref}/{chain_type}/{segment}.fasta"
    shell: "igseq tree before={input.before} after={input.after} {output}"

rule report_igdiscover_final_db:
    output: "analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.fasta"
    input: "analysis/igdiscover/{ref}/{chain_type}/{subject}/final/database/{segment}.fasta"
    shell: "cp {input} {output}"

### MINING-D

rule report_miningd_combo_tree:
    """Make a tree of all subjects' MINING-D output compared with all known D sequences."""
    output: "analysis/reporting/mining-d/all.nex"
    input: "analysis/reporting/mining-d/all.msa.fasta"
    shell: "igseq tree {input} {output}"

rule report_miningd_combo_msa:
    """Align all subjects' MINING-D output with all known D sequences."""
    output: "analysis/reporting/mining-d/all.msa.fasta"
    input:
        subject="analysis/reporting/mining-d/all.fasta",
        refs="analysis/reporting/mining-d/refs.fa"
    shell: "cat {input.refs} {input.subject} | igseq msa --input-format fa - {output}"

rule report_miningd_combo_fasta:
    """Combine all subjects' MINING-D outputs into one FASTA."""
    output: "analysis/reporting/mining-d/all.fasta"
    input: TARGET_REPORT_MININGD_FASTAS
    params:
        subjects = AVAILABLE_MININGD["subject"]
    run:
        with open(output[0], "wt") as f_out:
            for subject, path in zip(params.subjects, input):
                with open(path) as f_in:
                    for rec in SeqIO.parse(f_in, "fasta"):
                        rec.id = f"{subject}_{rec.id}"
                        SeqIO.write(rec, f_out, "fasta-2line")

rule report_miningd_tree:
    """Make a tree of each subject's MINING-D output compared with all known D sequences."""
    output: "analysis/reporting/mining-d/{subject}/{subject}.nex"
    input:
        msa="analysis/reporting/mining-d/{subject}/{subject}.msa.fasta",
        subject="analysis/reporting/mining-d/{subject}/{subject}.fasta",
        refs="analysis/reporting/mining-d/refs.fa"
    shell:
        """
            igseq tree \
                -C subject=#880000 -C refs=#000000 \
                -L subject=<(grep '^>' {input.subject} | cut -c 2- | cut -f 1 -d ' ') \
                -L refs=<(grep '^>' {input.refs} | cut -c 2- | cut -f 1 -d ' ') \
                {input.msa} {output}
        """

rule report_miningd_msa:
    """Align a subject's MINING-D output with all known D sequences."""
    output: "analysis/reporting/mining-d/{subject}/{subject}.msa.fasta"
    input:
        subject="analysis/reporting/mining-d/{subject}/{subject}.fasta",
        refs="analysis/reporting/mining-d/refs.fa"
    shell: "cat {input.refs} {input.subject} | igseq msa --input-format fa - {output}"

rule report_miningd_refs:
    """
    Make a combined D gene FASTA across rhesus refs.

    Output has one line per unqiue sequence across all references, with
    automated seq IDs, and IDs within individual references in the
    descriptions.
    """
    output: "analysis/reporting/mining-d/refs.fasta"
    run:
        from igseq.util import FILES, DATA
        from base64 import b32encode
        from hashlib import sha1
        refs = {
            "germ/rhesus/imgt/IGH/IGHD.fasta": "IMGT",
            "germ/rhesus/sonarramesh/IGH/IGHD.fasta": "Ramesh",
            "germ/rhesus/kimdb/IGH/IGHD.fasta": "KIMDB"}
        ref_seqs = defaultdict(list)
        for path in FILES:
            key = str(path.relative_to(DATA))
            if key not in refs:
                continue
            with open(path) as f_in:
                for rec in SeqIO.parse(f_in, "fasta"):
                    ref_seqs[str(rec.seq)].append((refs[key], rec.id))
        with open(output[0], "wt") as f_out:
            for seq, id_list in ref_seqs.items():
                # ID syntax inspired by 10.1016/j.immuno.2023.100025
                # except I'm making it deterministic (why make it random when you
                # can make it reproducible?)
                seqhash = sha1()
                seqhash.update(seq.encode("UTF8"))
                seqhash = b32encode(seqhash.digest()).decode("ASCII")
                seqid = "IGHD0-" + seqhash[:4] + "*00"
                seqdesc = " ".join([f"{ref}:{seqid}" for ref, seqid in id_list])
                f_out.write(f">{seqid} {seqdesc}\n{seq}\n")

rule report_miningd_output:
    output: "analysis/reporting/mining-d/{subject}/{subject}.fasta"
    input: "analysis/mining-d/{subject}.output.fasta"
    shell: "cp {input} {output}"
