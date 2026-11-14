#!/usr/bin/env python

"""
Summarize partis info per-lineage-group further, one row per lineage group.
"""

import re
import argparse
from collections import defaultdict
from csv import DictReader, DictWriter

DEFAULT_SORT_JUNCT_MIN=23
DEFAULT_SORT_TOTAL_MIN=2

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
        f"{prefix}d_call_source": "igblast" if chain == "heavy" else "",
        f"{prefix}cdr3_aa_length_min": max(0, (junction_aa_length_min or 0) - 2),
        f"{prefix}cdr3_aa_length_max": max(0, (junction_aa_length_max or 0) - 2),
        f"{prefix}junction_aa_v_min": juncts_min[0][1] if juncts_min else None,
        f"{prefix}junction_aa_v_max": juncts_max[0][1] if juncts_max else None,
        f"{prefix}junction_aa_length_min": junction_aa_length_min,
        f"{prefix}junction_aa_length_max": junction_aa_length_max,
        }
    if chain != "heavy":
        del attrs[f"{prefix}d_call"]
        del attrs[f"{prefix}d_call_source"]
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

def _prep_uca_cols(rows, chain):
    suffix = {"heavy": "", "light": "_light"}[chain]
    uca_ids = set()
    uca_seqs = set()
    diffs = (set(), set(), set())
    for row in rows:
        uca_ids.add(row[f"uca_sequence_id{suffix}"])
        uca_seqs.add(row[f"uca_sequence{suffix}"])
        tp = row["timepoint"]
        for idx in range(3):
            diffs_here = set(row[f"uca_cdr{idx+1}_aa_length_diff{suffix}"].split(" ")) - {""}
            diffs[idx].update({f"wk{tp}{diff}" for diff in diffs_here})
    diffs = tuple(" ".join(sorted(diff)) for diff in diffs)
    assert len(uca_ids) == 1
    assert len(uca_seqs) == 1
    uca_ids = uca_ids.pop()
    uca_seqs = uca_seqs.pop()
    out = {
        f"uca_sequence_id{suffix}": uca_ids,
        f"uca_sequence{suffix}": uca_seqs}
    out.update(
        {f"uca_cdr{idx+1}_aa_length_diff{suffix}": diff for idx, diff in enumerate(diffs)})
    return out

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

def _override_d_call(row_out, lin_attrs):
    # If we have a D call in the lineage attributes, use that instead of
    # whatever we gathered up from the post-partis summarizing via IgBLAST
    if (d_call := lin_attrs.get(row_out["lineage_group"], {}).get("DH")):
        row_out["d_call"] = d_call
        row_out["d_call_source"] = "lineage_attrs"

def _label_d_call_ref(row_out, d_refs):
    # there are some edge cases like where I've summarized differing IgBLAST
    # calls with "IGHD.../IGHD..." or where I've labeled VDDJ cases with
    # "IGHD...+IGHD..."
    d_call_refs = [d_refs.get(part) for part in re.split(r"[/\+]+", row_out["d_call"] or "")]
    out = d_call_refs[0] if len(set(d_call_refs)) == 1 else "/".join(d_call_refs)
    row_out["d_call_ref"] = out

def _prep_cols(info):
    # define categories, timepoints, output columns
    categories = sorted({row["category"] for row in info})
    timepoints = sorted({int(row["timepoint"]) for row in info})
    # [f"uca_cdr{idx+1}_aa_length_diff{suffix}"
    cols = [
        "lineage_group", "names",
        "uca_sequence_id", "uca_sequence"] + \
            [f"uca_cdr{idx+1}_aa_length_diff" for idx in range(3)] + [
        "uca_sequence_id_light", "uca_sequence_light"] + \
            [f"uca_cdr{idx+1}_aa_length_diff_light" for idx in range(3)] + [
        "lineage_attrs", "comments", "suspected_vddj"]
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
        "d_call_source",
        "d_call_ref",
        "cdr3_aa_length_min",
        "cdr3_aa_length_max",
        "junction_aa_v_max", # (ditto)
        "junction_aa_v_min", # (ditto)
        "junction_aa_length_min",
        "junction_aa_length_max"]
    cols += cols_heavy + [f"light_{col}" for col in cols_heavy if "d_call" not in col]
    return cols, categories

def _write(csv_out, cols, out):
    with open(csv_out, "w", encoding="ASCII") as f_out:
        writer = DictWriter(f_out, cols, lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)

def _load_d_refs(d_refs):
    # given list of refname=refpath pairs, set up dictionary mapping seq ID =>
    # ref name.  Requires all IDs to be unique across all refs.
    d_ref_for_id = {}
    if d_refs:
        for d_ref in d_refs:
            # parse custom ref name if present, otherwise just use path as name
            parts = d_ref.split("=", 1)
            if len(parts) == 1:
                parts = [parts[0], parts[0]]
            ref_name = parts[0]
            ref_path = parts[1]
            with open(ref_path, encoding="ASCII") as f_in:
                for line in f_in:
                    if not line.startswith(">"):
                        continue
                    seqid = line[1:].strip()
                    assert seqid not in d_ref_for_id, f"clashing ID for {seqid}"
                    d_ref_for_id[seqid] = ref_name
    return d_ref_for_id

def partis_lineages_summary(csv_in, csv_out, csv_in_lineage_attrs=None, *,
        paths_d_refs=None, sort_junct_min=DEFAULT_SORT_JUNCT_MIN, sort_total_min=DEFAULT_SORT_TOTAL_MIN):
    """Summarize partis info per-lineage-group further, one row per lineage group

    Sorting is first into those above or below a juction length cutoff, then
    those above or below a total lineage member count cutoff, then by
    decreasing junction length and member count numbers.

    sort_junct_min: at or above this heavy chain junction AA length (for max
                    observed across lineage) will be at top
    sort_total_min: next, at or above this for total member count will be at
                    top
    """
    with open(csv_in, encoding="ASCII") as f_in:
        info = list(DictReader(f_in))
    lin_attrs = {}
    if csv_in_lineage_attrs:
        with open(csv_in_lineage_attrs, encoding="ASCII") as f_in:
            lin_attrs = {row["Lineage"]: row for row in DictReader(f_in)}
    d_refs = _load_d_refs(paths_d_refs)
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
            "names": _condense_names(rows),
            "lineage_attrs": "Y" if group in lin_attrs else "",
            "comments": lin_attrs.get(group, {}).get("Comments"),
            "suspected_vddj": lin_attrs.get(group, {}).get("SuspectedVDDJ")}
        row_out.update(_prep_chain_attrs(rows, "heavy"))
        row_out.update(_prep_chain_attrs(rows, "light"))
        row_out.update(_prep_uca_cols(rows, "heavy"))
        row_out.update(_prep_uca_cols(rows, "light"))
        _add_summary_cols(row_out, rows, categories)
        _override_d_call(row_out, lin_attrs)
        _label_d_call_ref(row_out, d_refs)
        out.append(row_out)
    # sort with these things highest:
    #
    #  * >=whatever AA CDRH3
    #  * multiplets (that is, if sort_total_min=2)
    #  * longest junctions
    #  * highest member count
    out.sort(key=lambda row: (
        not row["junction_aa_length_max"] >= sort_junct_min,
        not row["total"] >= sort_total_min,
        -row["junction_aa_length_max"],
        -row["total"]))
    _write(csv_out, cols, out)

def main():
    """CLI for partis_lineages_summary"""
    parser = argparse.ArgumentParser()
    arg = parser.add_argument
    arg("input",
        help="CSV with per lineage group+timepoint+dataset category summary information")
    arg("output",
        help="CSV to write with one row per per lineage group")
    arg("-L", "--lineage-attrs", help="optional CSV with per-lineage metadata")
    arg("-D", "--d-ref", action="append", help="optional FASTA(s) for alternate D ref names")
    arg("-J", "--sort-junct-min", type=int, default=DEFAULT_SORT_JUNCT_MIN, help="")
    arg("-T", "--sort-total-min", type=int, default=DEFAULT_SORT_TOTAL_MIN, help="")
    args = parser.parse_args()
    partis_lineages_summary(
        args.input, args.output, args.lineage_attrs,
        paths_d_refs=args.d_ref,
        sort_junct_min=args.sort_junct_min,
        sort_total_min=args.sort_total_min)

if __name__ == "__main__":
    main()
