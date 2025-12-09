#!/usr/bin/env python

"""
Summarize partis info per-lineage-group further, one row per lineage group.
"""

import argparse
from collections import defaultdict
from csv import DictReader, DictWriter

def _condense_names(rows):
    # note isolate names if there aren't too many
    names = set()
    for row in rows:
        if "isolate" in row["category"]:
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

def get_junction_aa_lengths(rows, prefix, case):
    """Nonzero junction aa length mins or maxes across rows"""
    key = f"{prefix}junction_aa_length_{case}"
    return [int(row[key]) for row in rows if row[key]]

def get_juncts_for(rows, prefix, vcase):
    """Non-empty junctions paired with corresponding v identities for max/min case"""
    juncts = [(
        float(row[f"{prefix}v_identity_{vcase}"]) if row[f"{prefix}v_identity_{vcase}"] else 0,
        row[f"{prefix}junction_aa_v_{vcase}"]) for row in rows]
    juncts = [pair for pair in juncts if pair[1]]
    juncts.sort(reverse=vcase == "max")
    return juncts

def _prep_chain_attrs(rows, chain):
    prefix = {"heavy": "", "light": "light_"}[chain]
    # Minimum and maximum junction AA length across all sequences
    # across all categories and timepoints
    junction_aa_length_min = min(get_junction_aa_lengths(rows, prefix, "min"), default=None)
    junction_aa_length_max = max(get_junction_aa_lengths(rows, prefix, "max"), default=None)
    # The pair of junction AA sequences associated with the lowest V
    # identity and the highest V identity across all sequences across
    # all categories and timepoints
    juncts_min = get_juncts_for(rows, prefix, "min")
    juncts_max = get_juncts_for(rows, prefix, "max")
    attrs = {
        f"{prefix}v_family": _format_vj_family(rows, f"{prefix}v_family"),
        f"{prefix}j_family": _format_vj_family(rows, f"{prefix}j_family"),
        f"{prefix}v_identity_min": juncts_min[0][0] if juncts_min else None,
        f"{prefix}v_identity_max": juncts_max[0][0] if juncts_max else None,
        f"{prefix}d_call": _format_d_call(rows) if chain == "heavy" else "",
        f"{prefix}cdr3_aa_length_min": max(0, (junction_aa_length_min or 0) - 2),
        f"{prefix}cdr3_aa_length_max": max(0, (junction_aa_length_max or 0) - 2),
        f"{prefix}junction_aa_v_min": juncts_min[0][1] if juncts_min else None,
        f"{prefix}junction_aa_v_max": juncts_max[0][1] if juncts_max else None,
        f"{prefix}junction_aa_length_min": junction_aa_length_min,
        f"{prefix}junction_aa_length_max": junction_aa_length_max,
        }
    if chain != "heavy":
        del attrs[f"{prefix}d_call"]
    return attrs

def _format_vj_family(rows, key):
    # V or J family across everything (should only be one!  but if not, same
    # idea as for D, just no upper limit.)
    family = set()
    for row in rows:
        family = family | set(row[key].split("/"))
    family = family - {""}
    family = "/".join(sorted(family))
    return family

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
    cols_heavy = [
        "v_family",
        "j_family",
        "v_identity_max", # (first since max identity means min divergence)
        "v_identity_min", # (ditto)
        "d_call",
        "cdr3_aa_length_min",
        "cdr3_aa_length_max",
        "junction_aa_v_max", # (ditto)
        "junction_aa_v_min", # (ditto)
        "junction_aa_length_min",
        "junction_aa_length_max"]
    cols += cols_heavy + [f"light_{col}" for col in cols_heavy if col != "d_call"]
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
        row_out = {
            "lineage_group": group,
            "names": _condense_names(rows)}
        row_out.update(_prep_chain_attrs(rows, "heavy"))
        row_out.update(_prep_chain_attrs(rows, "light"))
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
