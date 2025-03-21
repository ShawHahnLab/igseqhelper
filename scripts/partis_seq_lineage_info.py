#!/usr/bin/env python

"""
Report sequences with partis clones overlapping with our isolates.
"""

import re
import gzip
import argparse
from collections import defaultdict
from csv import DictReader, DictWriter

def _load_isol_annots(csv_isol_annots):
    isolates = {}
    if csv_isol_annots:
        with open(csv_isol_annots, encoding="ASCII") as f_in:
            isolates = {row["Isolate"]: row for row in DictReader(f_in)}
    return isolates

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

def _prep_seq_lineage_info(airr_in, isolates, ngs_annots, igblast, keep_all):
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
                    if attrs.get("sequence") and \
                            annots.get("sequence") and \
                            attrs["sequence"] != annots["sequence"]:
                        raise ValueError(f"Sequence mismatch for {row['sequence_id']}")
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
    return out

def _assign_lineage_groups(out):
    clone_lineages = defaultdict(set) # clone ID -> set of lineages
    for row in out:
        clone_lineages[row["partis_clone_id"]].add(row["lineage"])
    # Clone IDs that include sequences without a lineage assigned will be used
    # to gather up *all* sequences referencing that clone ID, across whatever
    # lineages do happen to be assigned, into one lineage group.
    for row in out:
        lins = clone_lineages[row["partis_clone_id"]]
        row["lineage_group"] = ""
        if lins == {""}:
            # if no lineages were assigned to any of the sequences with this
            # clone ID assigned, just label it by the clone ID (if there is
            # one; really weird-looking sequences may not get one assigned)
            if row["partis_clone_id"]:
                row["lineage_group"] = "partis-" + row["partis_clone_id"]
        elif "" in lins:
            # otherwise, if it's a mix of assigned and blank lineages, use this
            # clone ID to group by all observed lineage names for this clone.
            # (Wait, should I also span across other clone IDs that overlap by
            # lineage name, too?  Not sure.)
            lineages = set(clone_lineages[row["partis_clone_id"]]) - {""}
            row["lineage_group"] = "/".join(sorted(lineages))
        else:
            # otherwise just use this row's one lineage as its group name,
            # ignoring partis' grouping
            row["lineage_group"] = row["lineage"]

def partis_seq_lineage_info(airr_in, csv_out, isol_annots=None,
        csv_ngs_annots=None, airr_in_igblast=None, *, keep_all=False):
    """Report sequences with partis clones overlapping with our isolates"""
    isolates = _load_isol_annots(isol_annots) # isolate names -> attrs
    ngs_annots = _load_ngs_annots(csv_ngs_annots) # seq ID -> custom attrs incl. Lineage
    igblast = _load_igblast_airr(airr_in_igblast) # seq ID -> IgBLAST attrs
    out = _prep_seq_lineage_info(airr_in, isolates, ngs_annots, igblast, keep_all)
    _assign_lineage_groups(out)
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
