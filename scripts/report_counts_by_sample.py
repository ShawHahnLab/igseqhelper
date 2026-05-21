#!/usr/bin/env python

"""
Make sample read count summary CSV (one row per sample)
"""

from pathlib import Path
from csv import DictReader, DictWriter
from argparse import ArgumentParser

ANALYSIS = Path("analysis")

def divide(val1, val2, fmt="{:.2f}"):
    """Divide val1 by val2 as floats and return formatted string."""
    try:
        num = float(val1)/float(val2)
    except (ValueError, TypeError, ZeroDivisionError):
        num = ""
    else:
        num = fmt.format(num)
    return num

def __report_counts_for(samp, runid, category):
    path = ANALYSIS/category/runid/f"{samp}.{category}.counts.csv"
    if path.exists():
        with open(path, encoding="ASCII") as f_in:
            reader = DictReader(f_in)
            for row in reader:
                if row["Item"] == "output":
                    return row["NumSeqs"]
    return ""

def __write(csv_out, rows_out, by_run):
    fieldnames = [
        "Run", "Subject", "Specimen", "CellType", "Type", "Sample",
        "CountsDemux", "CountsTrim", "CountsMerge", "CountsFilt",
        "CellCount", "RatioDemux", "RatioMerge", "RatioFilt"]
    with open(csv_out, "wt", encoding="ASCII") as f_out:
        writer = DictWriter(
            f_out,
            fieldnames=fieldnames,
            lineterminator="\n")
        writer.writeheader()
        # write per-run info
        for runid, rows in by_run["demux"].items():
            for row in rows:
                if row.get("Item") == "unassigned":
                    rows_out.append({
                        "Run": runid,
                        "Sample": "unassigned",
                        "CountsDemux": row["NumSeqs"]})
        for runid, rows in by_run["phix"].items():
            for row in rows:
                if row.get("Item") == "mapped":
                    rows_out.append({
                        "Run": runid,
                        "Sample": "unassigned.phix",
                        "CountsDemux": row["NumSeqs"]})
        rows_out = sorted(rows_out, key=lambda r: (str(r.get("Run")), str(r.get("Sample"))))
        writer.writerows(rows_out)

def report_counts_by_sample(samples_csv_in, specimens_csv_in, csv_out):
    """Make sample read count summary CSV (one row per sample)"""
    with open(samples_csv_in, encoding="ASCII") as f_in:
        samples = {row["Sample"]: row for row in DictReader(f_in)}
    with open(specimens_csv_in, encoding="ASCII") as f_in:
        specimens = {row["Specimen"]: row for row in DictReader(f_in)}
    for attrs in samples.values():
        attrs["SpecimenAttrs"] = specimens[attrs["Specimen"]]
    by_run = {"demux": {}, "phix": {}}
    rows_out = []
    for samp, attrs in samples.items():
        counts = {}
        path_demux = ANALYSIS/"demux"/attrs["Run"]/"demux.counts.csv"
        path_phix = ANALYSIS/"phix"/attrs["Run"]/"phix.counts.csv"
        counts["demux"] = ""
        if path_demux.exists():
            if attrs["Run"] not in by_run["demux"]:
                with open(path_demux, encoding="ASCII") as f_in:
                    rows = list(DictReader(f_in))
                by_run["demux"][attrs["Run"]] = rows
            for row in by_run["demux"][attrs["Run"]]:
                if row.get("Sample") == samp:
                    counts["demux"] = row["NumSeqs"]
        if path_phix.exists():
            if attrs["Run"] not in by_run["phix"]:
                with open(path_phix, encoding="ASCII") as f_in:
                    rows = list(DictReader(f_in))
                by_run["phix"][attrs["Run"]] = rows
        counts["trim"]  = __report_counts_for(samp, attrs["Run"], "trim")
        counts["merge"] = __report_counts_for(samp, attrs["Run"], "merge")
        counts["filt"]  = __report_counts_for(samp, attrs["Run"], "filt")
        rows_out.append({
            "Run": attrs["Run"],
            "Subject": attrs["SpecimenAttrs"]["Subject"],
            "Specimen": attrs["Specimen"],
            "CellType": attrs["SpecimenAttrs"]["CellType"],
            "Type": attrs["Type"],
            "Sample": samp,
            "CountsDemux": counts["demux"],
            "CountsTrim": counts["trim"],
            "CountsMerge": counts["merge"],
            "CountsFilt": counts["filt"],
            "CellCount": attrs["SpecimenAttrs"]["CellCount"],
            "RatioDemux": divide(counts["demux"], attrs["SpecimenAttrs"]["CellCount"]),
            "RatioMerge": divide(counts["merge"], attrs["SpecimenAttrs"]["CellCount"]),
            "RatioFilt": divide(counts["filt"], attrs["SpecimenAttrs"]["CellCount"])})
    __write(csv_out, rows_out, by_run)

def main():
    """CLI for report_counts_by_sample"""
    parser = ArgumentParser(description=__doc__)
    arg = parser.add_argument
    arg("samples", help="Path to samples metadata CSV")
    arg("specimens", help="Path to specimens metadata CSV")
    arg("output", help="Path for output CSV to write")
    args = parser.parse_args()
    report_counts_by_sample(args.samples, args.specimens, args.output)

if __name__ == "__main__":
    main()
