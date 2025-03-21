#!/usr/bin/env python

"""
Summarize partis info per-lineage-group, one row per group+timepoint+category.
"""

import re
import argparse
from collections import defaultdict
from csv import DictReader, DictWriter

COLS = [
    "lineage_group",
    "v_family",
    "v_identity_min",
    "v_identity_max",
    "d_call",
    "junction_aa_v_min",
    "junction_aa_v_max",
    "junction_aa_length_min",
    "junction_aa_length_max",
    "timepoint",
    "category",
    "total"]

def _load_seq_info(csv_in):
    with open(csv_in, encoding="ASCII") as f_in:
        seq_info = list(DictReader(f_in))
    groups = defaultdict(list)
    for row in seq_info:
        # group rows by lineage group, including those with none assigned
        # (in case all sequences were included in the input CSV)
        groups[row["lineage_group"]].append(row)
    return groups

def partis_lineages(csv_in, csv_out):
    """Summarize partis info per-lineage-group, one row per group+timepoint+category"""
    groups = _load_seq_info(csv_in)
    out = []
    for lineage_group, rows in groups.items():
        # within each lineage group, group rows by timepoint and category
        chunk = defaultdict(list)
        for row in rows:
            timepoint = int(row["timepoint"])
            category = row["category"]
            chunk[(timepoint, category)].append(row)
        for key, rows in chunk.items():
            v_family = ""
            d_call = ""
            junction_aa_length_min = min(int(row["junction_aa_length"]) for row in rows)
            junction_aa_length_max = max(int(row["junction_aa_length"]) for row in rows)
            juncts = None
            if lineage_group:
                v_family = "/".join(sorted({row["v_family"] for row in rows}))
                d_call = set()
                for row in rows:
                    d_call = d_call | set(re.split("[/,]", row["d_call"]))
                d_call = "/".join(sorted(d_call - {""}))
                junction_aa_length_min = min(int(row["junction_aa_length"]) for row in rows)
                junction_aa_length_max = max(int(row["junction_aa_length"]) for row in rows)
                juncts = [(
                    float(row["v_identity"]) if row["v_identity"] else 0,
                    row["junction_aa"]) for row in rows]
                juncts.sort()
            # set up output per group+timepoint+category
            out.append({
                "lineage_group": lineage_group,
                "v_family": v_family,
                "v_identity_min": f"{juncts[0][0]:.2f}" if juncts else "",
                "v_identity_max": f"{juncts[-1][0]:.2f}" if juncts else "",
                "d_call": d_call,
                "junction_aa_v_min": juncts[0][1] if juncts else "",
                "junction_aa_v_max": juncts[-1][1] if juncts else "",
                "junction_aa_length_min": junction_aa_length_min,
                "junction_aa_length_max": junction_aa_length_max,
                "timepoint": key[0],
                "category": key[1],
                "total": len(rows)})
    # total counts per lineage goup across timepoints and categories
    totals = defaultdict(int)
    for row in out:
        totals[row["lineage_group"]] += row["total"]
    # sort by lineage group by decreasing total count, then by increasing
    # timepoint, then by category
    def sorter(row):
        return (
            -totals[row["lineage_group"]],
            row["lineage_group"],
            row["timepoint"],
            row["category"])
    out.sort(key=sorter)
    with open(csv_out, "w", encoding="ASCII") as f_out:
        writer = DictWriter(f_out, COLS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)

def main():
    """CLI for partis_lineages"""
    parser = argparse.ArgumentParser()
    arg = parser.add_argument
    arg("input",
        help="CSV with per-sequence summary info derived from partis partitioning")
    arg("output",
        help="CSV to write with per lineage group+timepoint+dataset cateogry summary information")
    args = parser.parse_args()
    partis_lineages(args.input, args.output)

if __name__ == "__main__":
    main()
