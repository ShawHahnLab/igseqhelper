#!/usr/bin/env python

"""
Specimen read and SONAR cluster summary CSV (one row per specimen per cell+chain type)
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

def report_counts_by_specimen(counts_by_sample_csv, output_csv):
    """ Specimen read and SONAR cluster summary CSV (one row per specimen per cell+chain type)"""
    spec_info = {}
    with open(counts_by_sample_csv, encoding="ASCII") as f_in:
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
                sonar_counts = ANALYSIS/"reporting/counts/sonar"/ \
                    f"{row['Subject']}.{row['Type']}.{row['Specimen']}.csv"
                if sonar_counts.exists():
                    with open(sonar_counts, encoding="ASCII") as f_in_sonar:
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
    with open(output_csv, "wt", encoding="ASCII") as f_out:
        writer = DictWriter(
            f_out,
            fieldnames=fieldnames,
            lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def main():
    """CLI for report_counts_by_specimen"""
    parser = ArgumentParser(description=__doc__)
    arg = parser.add_argument
    arg("samples", help="Path to samples counts CSV")
    arg("output", help="Path for output CSV to write")
    args = parser.parse_args()
    report_counts_by_specimen(args.samples, args.output)

if __name__ == "__main__":
    main()
