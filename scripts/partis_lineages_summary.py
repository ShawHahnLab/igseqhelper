#!/usr/bin/env python

"""
Summarize partis info per-lineage-group further, one row per lineage group.
"""

import argparse
from collections import defaultdict
from csv import DictReader, DictWriter

def partis_lineages_summary(csv_in, csv_out):
    """Summarize partis info per-lineage-group further, one row per lineage group"""
    with open(csv_in, encoding="ASCII") as f_in:
        info = list(DictReader(f_in))
    # define categories, timepoints, output columns
    categories = sorted({row["category"] for row in info})
    timepoints = sorted({int(row["timepoint"]) for row in info})
    cols = ["lineage_group"]
    for category in categories:
        for timepoint in timepoints:
            cols.append(f"total_{category}_wk{timepoint}")
        cols.append(f"total_{category}")
        cols.append(f"earliest_{category}")
    cols.append("earliest")
    cols += [
        "v_family",
        "v_identity_min",
        "v_identity_max",
        "d_call",
        "junction_aa_v_min",
        "junction_aa_v_max",
        "junction_aa_length_min",
        "junction_aa_length_max"]
    # group by lineage group
    groups = defaultdict(list)
    for row in info:
        groups[row["lineage_group"]].append(row)
    # define output rows
    out = []
    for group, rows in groups.items():
        # D call(s) across everything (parsing and re-formatting any with
        # more than one to properly handle ["X/Y", "X"] -> "X/Y").  If
        # there are more than ten, throw up our hands and just put "???"
        d_call = set()
        for row in rows:
            d_call = d_call | set(row["d_call"].split("/"))
        d_call = d_call - {""}
        if len(d_call) > 10:
            d_call = "???"
        else:
            d_call = "/".join(sorted(d_call))
        # V family across everything (should only be one!  but if not, same
        # idea as for D, just no upper limit.)
        v_family = set()
        for row in rows:
            v_family = v_family | set(row["v_family"].split("/"))
        v_family = v_family - {""}
        v_family = "/".join(sorted(v_family))
        # Minimum and maximum junction AA length across all sequences
        # across all categories and timepoints
        junction_aa_length_min = min(int(row["junction_aa_length_min"]) for row in rows)
        junction_aa_length_max = max(int(row["junction_aa_length_max"]) for row in rows)
        # The pair of junction AA sequences associated with the lowest V
        # identity and the highest V identity across all sequences across
        # all categories and timepoints
        juncts_min = [(
            float(row["v_identity_min"]) if row["v_identity_min"] else 0,
            row["junction_aa_v_min"]) for row in rows]
        juncts_min.sort()
        juncts_max= [(
            float(row["v_identity_max"]) if row["v_identity_max"] else 0,
            row["junction_aa_v_max"]) for row in rows]
        juncts_max.sort(reverse=True)
        row_out = {
            "lineage_group": group,
            "v_family": v_family,
            "v_identity_min": juncts_min[0][0],
            "v_identity_max": juncts_max[0][0],
            "d_call": d_call,
            "junction_aa_v_min": juncts_min[0][1],
            "junction_aa_v_max": juncts_max[0][1],
            "junction_aa_length_min": junction_aa_length_min,
            "junction_aa_length_max": junction_aa_length_max,
            }
        totals = defaultdict(int)
        times = defaultdict(list)
        for row in rows:
            category = row["category"]
            timepoint = int(row["timepoint"])
            row_out[f"total_{category}_wk{timepoint}"] = int(row["total"])
            totals[category] += int(row["total"])
            times[category].append(timepoint)
        for category in categories:
            row_out[f"total_{category}"] = totals[category]
            row_out[f"earliest_{category}"] = min(times[category]) if times[category] else None
        mins = [min(vals) for vals in times.values() if vals]
        row_out["earliest"] = min(mins) if mins else None
        out.append(row_out)
    with open(csv_out, "w", encoding="ASCII") as f_out:
        writer = DictWriter(f_out, cols, lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)

def main():
    """CLI for partis_lineages_summary"""
    parser = argparse.ArgumentParser()
    arg = parser.add_argument
    arg("input",
        help="CSV with per lineage group+timepoint+dataset category summary information")
    arg("output",
        help="CSV to write with one row per per lineage group")
    args = parser.parse_args()
    partis_lineages_summary(args.input, args.output)

if __name__ == "__main__":
    main()
