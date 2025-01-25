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
        "CountsDemux", "CountsTrim", "CountsMerge", "CountsFilt",
        "CellCount", "RatioDemux", "RatioMerge", "RatioFilt"]
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
            samp_filt = report_counts_for(samp, attrs["Run"], "filt")
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
                "CountsFilt": samp_filt,
                "CellCount": attrs["SpecimenAttrs"]["CellCount"],
                "RatioDemux": divide(samp_demux, attrs["SpecimenAttrs"]["CellCount"]),
                "RatioMerge": divide(samp_merge, attrs["SpecimenAttrs"]["CellCount"]),
                "RatioFilt": divide(samp_filt, attrs["SpecimenAttrs"]["CellCount"])
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
                        "DemuxSeqs": [], "TrimSeqs": [], "MergeSeqs": [], "FiltSeqs": []}
                if row["CountsDemux"]:
                    spec_info[key]["DemuxSeqs"].append(int(row["CountsDemux"]))
                if row["CountsTrim"]:
                    spec_info[key]["TrimSeqs"].append(int(row["CountsTrim"]))
                if row["CountsMerge"]:
                    spec_info[key]["MergeSeqs"].append(int(row["CountsMerge"]))
                if row["CountsFilt"]:
                    spec_info[key]["FiltSeqs"].append(int(row["CountsFilt"]))
                # These take a while to crunch through so I'm doing that in a
                # separate rule, but still not explicitly giving it as input so
                # this rule will only use whatever's already on disk
                sonar_counts = Path(f"analysis/reporting/counts/sonar/{row['Subject']}.{row['Type']}.{row['Specimen']}.csv")
                if sonar_counts.exists():
                    with open(sonar_counts) as f_in_sonar:
                        reader = DictReader(f_in_sonar)
                        spec_info[key].update(next(reader))
    fieldnames = ["Subject", "Specimen", "CellType", "Type",
        "DemuxSeqs", "TrimSeqs", "MergeSeqs", "FiltSeqs",
        "CellCount", "RatioDemux", "RatioMerge", "RatioFilt",
        "SONARReads", "SONARGoodReads", "SONARClusteredReads", "SONARClusteredUnique",
        "LineageMembers"]
    rows = []
    for row in spec_info.values():
        row["DemuxSeqs"] = sum(row["DemuxSeqs"]) if row["DemuxSeqs"] else ""
        row["TrimSeqs"] = sum(row["TrimSeqs"]) if row["TrimSeqs"] else ""
        row["MergeSeqs"] = sum(row["MergeSeqs"]) if row["MergeSeqs"] else ""
        row["FiltSeqs"] = sum(row["FiltSeqs"]) if row["FiltSeqs"] else ""
        row["RatioDemux"] = divide(row["DemuxSeqs"], row["CellCount"])
        row["RatioMerge"] = divide(row["MergeSeqs"], row["CellCount"])
        row["RatioFilt"] = divide(row["FiltSeqs"], row["CellCount"])
        rows.append(row)
    rows = sorted(rows, key=lambda row: [row.get(field, "") for field in fieldnames])
    with open(output_csv, "wt") as f_out:
        writer = DictWriter(
            f_out,
            fieldnames=fieldnames,
            lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
