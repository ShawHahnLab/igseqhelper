#!/usr/bin/env python

"""
Run read count summary CSV (one row per run)
"""

import re
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

def __setup_run_row(run_info, row, quals):
    if row["Run"] not in run_info:
        runrow = {
            "Run": row["Run"],
            "UnassignedSeqs": "",
            "PhixSeqs": "",
            "SampleSeqs": []}
        if quals:
            for key in ("R1Q", "I1Q", "R2Q"):
                runrow[key] = quals.get(row["Run"], {}).get(key)
        run_info[row["Run"]] = runrow

def report_counts_by_run(counts_by_sample_csv, output_csv, quals_csv=None):
    """Run read count summary CSV (one row per run)"""
    run_info = {}
    with open(counts_by_sample_csv, encoding="ASCII") as f_in:
        rows = [row for row in DictReader(f_in) if row["Run"]]
    quals = {}
    if quals_csv:
        with open(quals_csv, encoding="ASCII") as f_in:
            quals = {row["Run"]: row for row in DictReader(f_in)}
    for row in rows:
        __setup_run_row(run_info, row, quals)
        path_run_counts = ANALYSIS/"reads"/row["Run"]/"getreads.counts.csv"
        if path_run_counts.exists():
            with open(path_run_counts, encoding="ASCII") as f_in:
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
    def run_date(attrs):
        match = re.match("(20)?([0-9]{6})[-_].*", attrs["Run"])
        txt1, txt2 = match.groups()
        txt1 = txt1 or "20"
        date = (txt1 + txt2[:2]) + "-" + txt2[2:4] + "-" + txt2[4:6]
        return date
    run_info = sorted(run_info.values(), key=run_date)
    with open(output_csv, "wt", encoding="ASCII") as f_out:
        writer = DictWriter(
            f_out,
            fieldnames=[
                "Run", "RawReads", "PassingFilter", "UnassignedSeqs", "PhixSeqs", "SampleSeqs",
                "UnassignedFraction", "PhixFraction", "R1Q", "I1Q", "R2Q"],
            lineterminator="\n")
        writer.writeheader()
        for row in run_info:
            row["SampleSeqs"] = sum(row["SampleSeqs"]) if row["SampleSeqs"] else ""
            parts = [row["UnassignedSeqs"], row["SampleSeqs"]]
            parts = [p for p in parts if p]
            row["TotalSeqs"] = sum(parts) if parts else ""
            row["UnassignedFraction"] = divide(row["UnassignedSeqs"], row["TotalSeqs"])
            row["PhixFraction"] = divide(row["PhixSeqs"], row["TotalSeqs"])
            del row["TotalSeqs"]
            writer.writerow(row)

def main():
    """CLI for report_counts_by_run"""
    parser = ArgumentParser(description=__doc__)
    arg = parser.add_argument
    arg("samples", help="Path to samples counts CSV")
    arg("--qualsummary", help="Path to optional per-run quality metrics CSV")
    arg("output", help="Path for output CSV to write")
    args = parser.parse_args()
    report_counts_by_run(args.samples, args.output, args.qualsummary)

if __name__ == "__main__":
    main()
