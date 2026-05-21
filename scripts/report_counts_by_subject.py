#!/usr/bin/env python

"""
Subject read and SONAR cluster summary CSV (one row per subject per cell+chain type)
"""

from csv import DictReader, DictWriter
from collections import defaultdict
from argparse import ArgumentParser

def report_counts_by_subject(counts_by_specimen_csv, output_csv):
    """Condense counts by specimen (by cell type+chain) down to by subject (by cell type+chain)"""
    # some input columns will be used for grouping, and some will be used to
    # total across
    keykeys = ("Subject", "CellType", "Type")
    totalkeys = (
        "DemuxSeqs", "TrimSeqs", "MergeSeqs", "FiltSeqs", "CellCount",
        "SONARReads", "SONARGoodReads", "SONARClusteredReads", "SONARClusteredUnique",
        "LineageMembers")
    # group input rows based on the columns given above
    grouped = defaultdict(list)
    with open(counts_by_specimen_csv, encoding="ASCII") as f_in:
        for row in DictReader(f_in):
            key = tuple(row[k] for k in keykeys)
            grouped[key].append(row)
    # Flatten into one row per group
    out = []
    for key, group in grouped.items():
        def total(attr, group):
            if all(row[attr] in (None, "") for row in group):
                return None
            return sum(int(row[attr] if row[attr] else 0) for row in group)
        # Put the grouping keys back to start with, followed by the number of
        # specimens per group, followed by the totals across specimens
        row_out = dict(zip(keykeys, key))
        row_out["Specimens"] = len(group)
        row_out.update({key: total(key, group) for key in totalkeys})
        out.append(row_out)
    with open(output_csv, "w", encoding="ASCII") as f_out:
        writer = DictWriter(f_out, out[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)

def main():
    """CLI for report_counts_by_subject"""
    parser = ArgumentParser(description=__doc__)
    arg = parser.add_argument
    arg("specimens", help="Path to specimens counts CSV")
    arg("output", help="Path for output CSV to write")
    args = parser.parse_args()
    report_counts_by_subject(args.specimens, args.output)

if __name__ == "__main__":
    main()
