from csv import DictWriter, DictReader

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
    rows_out = []
    with open(csv_out, "wt") as f_out:
        writer = DictWriter(
            f_out,
            fieldnames=fieldnames,
            lineterminator="\n")
        writer.writeheader()
        for samp, attrs in SAMPLES.items():
            path_demux = Path("analysis/demux")/attrs["Run"]/"demux.counts.csv"
            samp_demux = ""
            if path_demux.exists():
                if attrs["Run"] not in demux:
                    with open(path_demux) as f_in:
                        rows = list(DictReader(f_in))
                    demux[attrs["Run"]] = rows
                for row in demux[attrs["Run"]]:
                    if row.get("Sample") == samp:
                        samp_demux = row["NumSeqs"]
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
                        "SampleSeqs": []}
                if row["Sample"] == "unassigned":
                    if row["CountsDemux"]:
                        run_info[row["Run"]]["UnassignedSeqs"] = int(row["CountsDemux"])
                elif row["Sample"] == "unassigned.phix":
                    pass
                else:
                    if row["CountsDemux"]:
                        run_info[row["Run"]]["SampleSeqs"].append(int(row["CountsDemux"]))
    with open(output_csv, "wt") as f_out:
        writer = DictWriter(f_out, fieldnames=["Run", "UnassignedSeqs", "SampleSeqs", "TotalSeqs", "Ratio"], lineterminator="\n")
        writer.writeheader()
        for row in run_info.values():
            row["SampleSeqs"] = sum(row["SampleSeqs"]) if row["SampleSeqs"] else ""
            parts = [row["UnassignedSeqs"], row["SampleSeqs"]]
            parts = [p for p in parts if p]
            if parts:
                row["TotalSeqs"] = sum(parts)
            row["Ratio"] = divide(row["UnassignedSeqs"], row["SampleSeqs"])
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
