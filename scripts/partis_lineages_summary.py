#!/usr/bin/env python

"""
Summarize partis info per-lineage-group further, one row per lineage group.
"""

import argparse
from collections import defaultdict
from csv import DictReader, DictWriter

def _condense_names(rows):
    names = set()
    for row in rows:
        if row["category"] == "isolate":
            if row["names"] == "(many)":
                names = "(many)"
                break
            names |= set(row["names"].split("/"))
    else:
        names = sorted(names)
        # TODO would be nice to handle things like
        # CC01-40/CC01-41/.../CC01-45/CC01-46
        # as
        # "CC01-40 through CC01-46"
        names = "/".join(names) if len(names) < 10 else "(many)"
    return names

def _prep_row_out(group, rows):
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
    # note isolate names if there aren't too many
    names = _condense_names(rows)
    row_out = {
        "lineage_group": group,
        "names": names,
        "v_family": _format_v_family(rows),
        "v_identity_min": juncts_min[0][0],
        "v_identity_max": juncts_max[0][0],
        "d_call": _format_d_call(rows),
        "junction_aa_v_min": juncts_min[0][1],
        "junction_aa_v_max": juncts_max[0][1],
        "junction_aa_length_min": junction_aa_length_min,
        "junction_aa_length_max": junction_aa_length_max,
        }
    return row_out

def _format_v_family(rows):
    # V family across everything (should only be one!  but if not, same
    # idea as for D, just no upper limit.)
    v_family = set()
    for row in rows:
        v_family = v_family | set(row["v_family"].split("/"))
    v_family = v_family - {""}
    v_family = "/".join(sorted(v_family))
    return v_family

def _format_d_call(rows):
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
    return d_call

def _add_summary_cols(row_out, rows, categories):
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
        row_out[f"latest_{category}"] = max(times[category]) if times[category] else None
    row_out["total"] = sum(totals.values())
    mins = [min(vals) for vals in times.values() if vals]
    row_out["earliest"] = min(mins) if mins else None
    maxes = [max(vals) for vals in times.values() if vals]
    row_out["latest"] = max(maxes) if maxes else None

def _prep_cols(info):
    # define categories, timepoints, output columns
    categories = sorted({row["category"] for row in info})
    timepoints = sorted({int(row["timepoint"]) for row in info})
    cols = ["lineage_group", "names"]
    for category in categories:
        for timepoint in timepoints:
            # member total at this timepoint for this category
            cols.append(f"total_{category}_wk{timepoint}")
        # member total and earliest and latest timepoint of this category
        cols.append(f"total_{category}")
        cols.append(f"earliest_{category}")
        cols.append(f"latest_{category}")
    # member total and earliest and latest timepoint across all categories
    cols.append("total")
    cols.append("earliest")
    cols.append("latest")
    cols += [
        "v_family",
        "v_identity_min",
        "v_identity_max",
        "d_call",
        "junction_aa_v_min",
        "junction_aa_v_max",
        "junction_aa_length_min",
        "junction_aa_length_max"]
    return cols, categories

def partis_lineages_summary(csv_in, csv_out):
    """Summarize partis info per-lineage-group further, one row per lineage group"""
    with open(csv_in, encoding="ASCII") as f_in:
        info = list(DictReader(f_in))
    cols, categories = _prep_cols(info)
    # group by lineage group
    groups = defaultdict(list)
    for row in info:
        groups[row["lineage_group"]].append(row)
    # define output rows
    out = []
    for group, rows in groups.items():
        row_out = _prep_row_out(group, rows)
        _add_summary_cols(row_out, rows, categories)
        out.append(row_out)
    # sort with these things highest:
    #
    #  * >=20 AA CDRH3
    #  * multiplets
    #  * longest junctions
    #  * highest member count
    out.sort(key=lambda row: (
        not row["junction_aa_length_max"] >= 22,
        not row["total"] > 1,
        -row["junction_aa_length_max"],
        -row["total"]))
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
