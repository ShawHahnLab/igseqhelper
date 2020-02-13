"""
Summarizing and reporting helper functions.
"""
import csv
from igseq.data import load_csv

def counts_sample_summary(
        file_counts_in, counts_sample_summary_out, samples, runs, samples_per_run):
    """Take a list of sequence counts for fastq.gz files and summarize per-sample."""
    counts = load_csv(file_counts_in)
    rows = []
    for runid in runs.keys():
        for sample in samples_per_run[runid] + ["unassigned"]:
            key = "counts/demux/{run}/{sample}.I1.fastq.gz.counts".format(
                run=runid,
                sample=sample)
            cts = int(counts[key]["NumSequences"])
            try:
                cellcount = samples[sample]["SpecimenAttrs"]["CellCount"]
                ratio = divide(cts, cellcount)
            except KeyError:
                cellcount = ""
                ratio = ""
            rows.append({
                "Run": runid,
                "Sample": sample,
                "CellCount": cellcount,
                "NumSequences": cts,
                "Ratio": ratio})
    with open(counts_sample_summary_out, "wt") as f_out:
        writer = csv.DictWriter(
            f_out,
            fieldnames=["Run", "Sample", "NumSequences", "CellCount", "Ratio"])
        writer.writeheader()
        writer.writerows(rows)

def counts_run_summary(counts_sample_summary_in, counts_run_summary_out):
    """Take a per-sample summary of sequence counts and make a per-run version."""
    with open(counts_sample_summary_in) as f_in:
        reader = csv.DictReader(f_in)
        counts_by_sample = list(reader)
    counts_by_run = {}
    for row_in in counts_by_sample:
        runid = row_in["Run"]
        try:
            cts = counts_by_run[runid]
        except KeyError:
            counts_by_run[runid] = {"samples": 0, "unassigned": 0}
            cts = counts_by_run[runid]
        if row_in["Sample"] == "unassigned":
            key = "unassigned"
        else:
            key = "samples"
        cts[key] += int(row_in["NumSequences"])
    rows = []
    for run_name, run_attrs in counts_by_run.items():
        try:
            ratio = divide(run_attrs["unassigned"], run_attrs["samples"])
        except KeyError:
            ratio = ""
        rows.append([
            run_name,
            run_attrs.get("unassigned", ""),
            run_attrs.get("samples", ""),
            run_attrs.get("unassigned", "") + run_attrs.get("samples", ""),
            ratio])
    with open(counts_run_summary_out, "wt") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["Run", "UnassignedSeqs", "SampleSeqs", "TotalSeqs", "Ratio"])
        writer.writerows(rows)

def divide(val1, val2, fmt="{:.6f}"):
    """Divide val1 by val2 as floats and return formatted string."""
    num = float(val1)/float(val2)
    num = fmt.format(num)
    return num
