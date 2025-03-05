#!/usr/bin/env python

"""
Report sequences with partis clones overlapping with our isolates.
"""

import re
import gzip
import argparse
from collections import defaultdict
from csv import DictReader, DictWriter

def _load_ngs_annots(csv_ngs_annots):
    ngs_annots = {}
    if csv_ngs_annots:
        with open(csv_ngs_annots, encoding="ASCII") as f_in:
            for row in DictReader(f_in):
                ngs_annots[row["sequence_id"]] = row
    return ngs_annots

def _load_igblast_airr(airr_path):
    if airr_path:
        with gzip.open(airr_path, "rt", encoding="ASCII") as f_in:
            return {row["sequence_id"]: row for row in DictReader(f_in, delimiter="\t")}
    return {}

def partis_seq_lineage_info(airr_in, csv_out, isol_annots,
        csv_ngs_annots=None, airr_in_igblast=None, *, keep_all=False):
    # isolate -> atributes
    with open(isol_annots, encoding="ASCII") as f_in:
        isolates = {row["Isolate"]: row for row in DictReader(f_in)}
    # if given, load NGS sequence ID -> attributes including Lineage name
    ngs_annots = _load_ngs_annots(csv_ngs_annots)
    # if given, load IgBLAST AIRR as dictionary by sequence ID
    igblast = _load_igblast_airr(airr_in_igblast)
    # clone ID -> AIRR rows (need all to decide what to keep later)
    clones = defaultdict(list)
    # clone IDs of interest for our isolates
    cloneids = set()
    with open(airr_in, encoding="ASCII") as f_in:
        for row in DictReader(f_in, delimiter="\t"):
            clones[row["clone_id"]].append(row)
            if row["sequence_id"] in isolates:
                cloneids.add(row["clone_id"])
    # include everything that's listed under any of those clone IDs of
    # interest.  Each sequence can have one clone ID from partis and one
    # (if it's in our isolate metadata) Lineage assigned from us.
    out = []
    for cloneid, rows in clones.items():
        if cloneid in cloneids or keep_all:
            # keep all for this clone, or everything if specified
            for row in rows:
                timepoint_seqid = re.match("wk([0-9]+)-.*", row["sequence_id"])
                if timepoint_seqid:
                    timepoint_seqid = timepoint_seqid.group(1)
                # if there are IgBLAST-provided annotations, use those, but if
                # not, use what's already in this row.  Note this will also use
                # the sequence from IgBLAST if available since partis seems to
                # pad it with N for some reason
                annots = igblast.get(row["sequence_id"], row)
                v_family = re.match("(IG[HKL]V[0-9]+)", annots["v_call"])
                v_family = v_family.group(1) if v_family else ""
                # Is sequence from isolates?
                isol_attrs = isolates.get(row["sequence_id"], {})
                category = "isolate" if isol_attrs else "ngs"
                timepoint = isol_attrs.get("Timepoint", timepoint_seqid)
                lineage = isol_attrs.get("Lineage")
                # Or, do we have NGS seq annotations for it?
                if row["sequence_id"] in ngs_annots:
                    attrs = ngs_annots[row["sequence_id"]]
                    # sanity check with sequence content if present
                    if attrs.get("sequence") and attrs["sequence"] != row["sequence"]:
                        raise ValueError
                    if not lineage:
                        lineage = attrs.get("Lineage")
                row_out = {
                    "sequence_id": row["sequence_id"],
                    "sequence": annots["sequence"],
                    "v_family": v_family,
                    "v_identity": annots["v_identity"],
                    "d_call": annots["d_call"],
                    "junction_aa": annots["junction_aa"],
                    "junction_aa_length": len(annots["junction_aa"]),
                    "timepoint": timepoint,
                    "category": category,
                    "partis_clone_id": row["clone_id"] or "",
                    "lineage": lineage or ""}
                out.append(row_out)

    # below was previously in separate rule


    lineage_clones = defaultdict(set) # lineage -> set of clone IDs
    clone_lineages = defaultdict(set) # clone ID -> set of lineages
    # what clones are associated with what lineage?
    for row in out:
        clone_lineages[row["partis_clone_id"]].add(row["lineage"])
        lineage_clones[row["lineage"]].add(row["partis_clone_id"])
        #if row["partis_clone_id"] and row["lineage"]:
        #    clone_lineages[row["partis_clone_id"]].add(row["lineage"])
        #    lineage_clones[row["lineage"]].add(row["partis_clone_id"])
    for row in out:
        # group rows by lineage involved for that clone ID.  Hopefully just
        # one!  But for any edge cases with more than one, we'll note the
        # combo as a single string.  For rows that don't have the lineage
        # defined, leave that empty string out.
        lineages = ({row["lineage"]} | set(clone_lineages[row["partis_clone_id"]])) - {""}
        #lineages = {lineage for lineage in lineages if lineage}
        row["lineage_group"] = "/".join(sorted(lineages))
    out.sort(key = lambda row: (
        row["lineage_group"] == "",
        row["partis_clone_id"] == "",
        row["lineage_group"],
        row["lineage"],
        row["partis_clone_id"],
        row["sequence_id"]))
    keys = [
        "lineage_group",
        "sequence_id",
        "sequence",
        "v_family",
        "v_identity",
        "d_call",
        "junction_aa",
        "junction_aa_length",
        "timepoint",
        "category",
        "partis_clone_id",
        "lineage"]
    with open(csv_out, "w", encoding="ASCII") as f_out:
        writer = DictWriter(f_out, keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)

def main():
    """CLI for seq_lineage_info"""
    parser = argparse.ArgumentParser()
    arg = parser.add_argument
    arg("input", help="Partis AIRR TSV with clone_id")
    arg("output", help="CSV to write with summary sequence and lineage information")
    arg("-i", "--isolates", help="CSV with Isolate info")
    arg("-n", "--ngs-annotations", help="optional CSV with Lineage info for known NGS sequences")
    arg("-A", "--igblast-airr", help="optional AIRR tsv.gz from IgBLAST to prefer for annotations")
    arg("-a", "--all", action="store_true",
        help="keep all sequences or only those belonging to clones that also include isolates?")
    args = parser.parse_args()
    partis_seq_lineage_info(
        args.input, args.output, args.isolates,
        args.ngs_annotations, args.igblast_airr, keep_all=args.all)

if __name__ == "__main__":
    main()
