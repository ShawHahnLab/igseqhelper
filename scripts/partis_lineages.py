#!/usr/bin/env python

"""
Summarize partis info per-lineage-group, one row per group+timepoint+category.
"""

import re
import argparse
from collections import defaultdict
from csv import DictReader, DictWriter

COLS_HEAVY =  [
    "v_family",
    "j_family",
    "v_identity_min",
    "v_identity_max",
    "d_call",
    "junction_aa_v_min",
    "junction_aa_v_max",
    "junction_aa_length_min",
    "junction_aa_length_max"]
COLS_LIGHT = [f"light_{col}" for col in COLS_HEAVY if col != "d_call"]
COLS = [
    "lineage_group",
    "names"] + COLS_HEAVY + COLS_LIGHT + [
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

def condense_families(rows, key):
    return "/".join(sorted({row[key] for row in rows if row[key]}))

def get_junction_aa_lengths(rows, prefix):
    """Nonzero junction aa lengths"""
    return [r[f"{prefix}junction_aa_length"] for r in rows if r[f"{prefix}junction_aa_length"]]

def get_juncts(rows, prefix):
    """Non-empty junctions paired with corresponding v identities"""
    juncts = [(
        float(row[f"{prefix}v_identity"]) if row[f"{prefix}v_identity"] else 0,
        row[f"{prefix}junction_aa"]) for row in rows]
    juncts = [pair for pair in juncts if pair[1]]
    juncts.sort()
    return juncts

def _get_chain_attrs(rows, chain, lineage_group):
    prefix = {"heavy": "", "light": "light_"}[chain]
    v_family = ""
    j_family = ""
    d_call = ""
    junctlens = get_junction_aa_lengths(rows, prefix)
    junction_aa_length_min = min(junctlens, default=0)
    junction_aa_length_max = max(junctlens, default=0)
    juncts = None
    if lineage_group:
        v_family = condense_families(rows, f"{prefix}v_family")
        j_family = condense_families(rows, f"{prefix}j_family")
        d_call = set()
        if chain == "heavy":
            for row in rows:
                d_call = d_call | set(re.split("[/,]", row["d_call"]))
            d_call = "/".join(sorted(d_call - {""}))
        juncts = get_juncts(rows, prefix)
    attrs = {
        f"{prefix}v_family": v_family,
        f"{prefix}j_family": j_family,
        f"{prefix}v_identity_min": f"{juncts[0][0]:.2f}" if juncts else "",
        f"{prefix}v_identity_max": f"{juncts[-1][0]:.2f}" if juncts else "",
        f"{prefix}d_call": d_call, # only keep below if heavy
        f"{prefix}junction_aa_v_min": juncts[0][1] if juncts else "",
        f"{prefix}junction_aa_v_max": juncts[-1][1] if juncts else "",
        f"{prefix}junction_aa_length_min": junction_aa_length_min,
        f"{prefix}junction_aa_length_max": junction_aa_length_max}
    if chain != "heavy":
        del attrs[f"{prefix}d_call"]
    return attrs

def _group_rows(rows):
    chunk = defaultdict(list)
    for row in rows:
        timepoint = int(row["timepoint"])
        category = row["category"]
        chunk[(timepoint, category)].append(row)
    return chunk

def partis_lineages(csv_in, csv_out):
    """Summarize partis info per-lineage-group, one row per group+timepoint+category"""
    groups = _load_seq_info(csv_in)
    out = []
    for lineage_group, rows in groups.items():
        # within each lineage group, group rows by timepoint and category
        for key, rows in _group_rows(rows).items():
            # set up output per group+timepoint+category
            # note original IDs if there aren't too many of them
            names = ""
            names = {row["sequence_id_original"] for row in rows}
            names = "/".join(sorted(names)) if len(names) < 10 else "(many)"
            row_out = {
                "lineage_group": lineage_group,
                "names": names}
            row_out.update(_get_chain_attrs(rows, "heavy", lineage_group))
            row_out.update(_get_chain_attrs(rows, "light", lineage_group))
            row_out.update({
                "timepoint": key[0],
                "category": key[1],
                "total": len(rows)})
            out.append(row_out)
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
