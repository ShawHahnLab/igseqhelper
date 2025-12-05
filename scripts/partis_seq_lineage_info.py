#!/usr/bin/env python

"""
Report sequence and lineage info from partis, merging metadata for our seqs+isolates.
"""

import re
import gzip
import argparse
from collections import defaultdict
from csv import DictReader, DictWriter

def _load_metadata(csv_path, key=None):
    things = {}
    if csv_path:
        first = lambda row: list(row.keys())[0]
        with open(csv_path, encoding="ASCII") as f_in:
            things = {row[key or first(row)]: row for row in DictReader(f_in)}
    return things

def _load_igblast_airr(airr_path):
    if airr_path:
        with gzip.open(airr_path, "rt", encoding="ASCII") as f_in:
            return {row["sequence_id"]: row for row in DictReader(f_in, delimiter="\t")}
    return {}

def _load_clones_from_partis_airr(airr_in, metadata, keep_all):
    # clone ID -> AIRR rows (need all to decide what to keep later)
    clones = defaultdict(list)
    # clone IDs of interest for our isolates
    cloneids = set()
    with open(airr_in, encoding="ASCII") as f_in:
        for row in DictReader(f_in, delimiter="\t"):
            clones[row["clone_id"]].append(row)
            if not keep_all:
                if row["sequence_id"] in metadata["isolates"]:
                    cloneids.add(row["clone_id"])
    if keep_all:
        cloneids = None
    return clones, cloneids

def __infer_basics_from_metadata(seqid_in, metadata):
    category, seqid = seqid_in.split("-", 1)
    item = ""
    # isolate metadata is per isolate, the others are per item
    if category == "isolate":
        attrs = metadata.get(f"{category}s", {}).get(seqid, {})
    else:
        item, seqid = seqid.split("-", 1)
        attrs = metadata.get(f"{category}s", {}).get(item, {})
    timepoint_seqid = re.match("wk([0-9]+)-.*", seqid_in)
    if timepoint_seqid:
        timepoint_seqid = timepoint_seqid.group(1)
    row_out = {
        "sequence_id": seqid_in,
        "sequence_id_original": seqid,
        "category": category,
        "item": item,
        "timepoint": attrs.get("Timepoint", timepoint_seqid),
        "notes": "",
        "lineage": ""}
    # Catch the special case of seqset entries that have been added to the
    # isolates table
    if category == "seqset":
        isolate_alt_name = attrs["Subject"] + "-wk" + attrs["Timepoint"] + "-" + seqid
        isolate_map = {r["AltName"]: r for r in metadata["isolates"].values() if r["AltName"]}
        if isolate_alt_name in isolate_map:
            print(f"TODO handle seqset entry as isolate {isolate_alt_name}")
            # So, should I just overwrite the seqsets info with isolates info?
            # Maybe?  Could run into duplicates though
    if attrs.get("Skip") == "TRUE":
        # skip entries if that's noted in their metadata (isolates
        # we don't want to actually analyze, basically); generally
        # they shouldn't get this far anyway, but if so, we'll
        # exclude them now
        return row_out
    # A bit of post-processing on the category labels
    # (using the shorthand "ngs" for per-specimen material, and
    # including a more specific suffix for cases where a
    # particular preparation method was noted, like 10x)
    if row_out["category"] == "specimen":
        row_out["category"] = "ngs"
    if attrs.get("Method"):
        # (oh except don't both with a suffix for Sanger isolates;
        # that can just be left implicit, since most are Sanger)
        if not (row_out["category"] == "isolate" and attrs["Method"] == "Sanger"):
            row_out["category"] += "_" + attrs["Method"]
    # (Isolate lineage names that include the keyword "unassigned"
    # in the name will be interpreted as placeholders and ignored.
    # Other categories won't have a Lineage explicitly provided
    # anyway.)
    row_out["lineage"] = attrs.get("Lineage", "")
    if "unassigned" in row_out["lineage"]:
        row_out["lineage"] = ""
    return row_out

def __include_igblast_attrs(row_out, igblast):
    # if there are IgBLAST-provided annotations, use those, but if not, use
    # what's already in this row.  Note this will also use the sequence from
    # IgBLAST if available since partis seems to pad it with N for some reason
    igblast_attrs = igblast.get(row_out["sequence_id"], row_out)
    v_family = re.match("(IG[HKL]V[0-9]+)", igblast_attrs["v_call"])
    v_family = v_family.group(1) if v_family else ""
    row_out.update({
        "sequence": igblast_attrs["sequence"],
        "v_family": v_family,
        "v_identity": igblast_attrs["v_identity"],
        "d_call": igblast_attrs["d_call"],
        "junction_aa": igblast_attrs["junction_aa"],
        "junction_aa_length": len(igblast_attrs["junction_aa"])})

def __include_ngs_attrs(row_out, ngs_annots):
    # Do we have NGS seq annotations for it?  If so, take the lineage from
    # there if not otherwise specified
    if row_out["category"] == "ngs":
        ngsid = "wk" + row_out["timepoint"] + "-" + row_out["sequence_id_original"]
        if ngsid in ngs_annots:
            ngs_attrs = ngs_annots[ngsid]
            # sanity check with sequence content if present
            if ngs_attrs.get("sequence") and \
                    row_out.get("sequence") and \
                    ngs_attrs["sequence"] not in row_out["sequence"]:
                print(row_out["sequence"])
                print(ngs_attrs["sequence"])
                raise ValueError(f"Sequence mismatch for {row_out['sequence_id']}")
            if not row_out["lineage"]:
                row_out["lineage"] = ngs_attrs.get("Lineage", "")

def _prep_seq_lineage_info(clones, metadata, ngs_annots, igblast, cloneids):
    # include everything that's listed under any of those clone IDs of
    # interest, if defined.  Each sequence can have one clone ID from partis
    # and one (if it's in our isolate metadata) Lineage assigned from us.
    out = []
    for cloneid, rows in clones.items():
        if cloneids is None or cloneid in cloneids:
            # keep all for this clone, or everything if specified
            for row in rows:
                # start of by figuring out what sort of sequence this is, and
                # its metadata, from the sequence ID.  Category of None implies
                # skip this one entirely.
                row_out = __infer_basics_from_metadata(row["sequence_id"], metadata)
                if row_out["category"] is None:
                    continue
                # Add additional information with the help of IgBLAST output
                # and (if applicable) manually-defined info on NGS sequences
                __include_igblast_attrs(row_out, igblast)
                __include_ngs_attrs(row_out, ngs_annots)
                row_out.update({
                    "partis_clone_id": row["clone_id"] or "",
                    "lineage": row_out["lineage"] or ""})
                out.append(row_out)
    return out

def _assign_lineage_groups(out):
    # For each partis clone ID, note the set of all corresponding lineages we
    # have manually assigned from any data source
    clone_lineages = defaultdict(set)
    for row in out:
        clone_lineages[row["partis_clone_id"]].add(row["lineage"])
    # Clone IDs that include any sequences without a lineage assigned will be
    # used to gather up *all* sequences referencing that clone ID, across
    # whatever lineages do happen to be assigned, into one lineage group.
    for row in out:
        lins = clone_lineages[row["partis_clone_id"]]
        row["lineage_group_category"] = ""
        row["lineage_group"] = ""
        if lins == {""}:
            # If no lineages were assigned to any of the sequences with this
            # clone ID assigned, just label it by the clone ID (if there is
            # one; really weird-looking sequences may not get a clone ID
            # assigned at all)
            row["lineage_group_category"] = "none"
            if row["partis_clone_id"]:
                row["lineage_group"] = "partis-" + row["partis_clone_id"]
                row["lineage_group_category"] = "automatic"
        elif "" in lins:
            # Otherwise, if it's a mix of assigned and blank lineages, use this
            # clone ID to group by all observed lineage names for this clone.
            # (TODO Wait, should I also span across other clone IDs that
            # overlap by lineage name, too?  Yeah probably.  Currently this is
            # set up so that we could end with some things labeled "linA" and
            # others "linA/linB".  But good enough for now.)
            lineages = set(clone_lineages[row["partis_clone_id"]]) - {""}
            row["lineage_group"] = "/".join(sorted(lineages))
            row["lineage_group_category"] = "partis-grouped"
        else:
            # otherwise just use this row's one lineage as its group name,
            # ignoring partis' grouping
            row["lineage_group"] = row["lineage"]
            row["lineage_group_category"] = "manual"

def partis_seq_lineage_info(
        airr_in, csv_out,
        metadata_isolates=None, metadata_specimens=None, metadata_seqsets=None,
        csv_ngs_annots=None, airr_in_igblast=None, *, keep_all=False):
    """Report sequences with partis clones overlapping with our isolates"""
    # name -> attrs
    metadata = {
        "isolates": _load_metadata(metadata_isolates),
        "specimens": _load_metadata(metadata_specimens),
        "seqsets": _load_metadata(metadata_seqsets),
        }
    # seq ID here -> custom attrs incl. Lineage
    ngs_annots = _load_metadata(csv_ngs_annots, "sequence_id")
     # seq ID here -> IgBLAST attrs
    igblast_annots = _load_igblast_airr(airr_in_igblast)
    clones, cloneids = _load_clones_from_partis_airr(airr_in, metadata, keep_all)
    out = _prep_seq_lineage_info(clones, metadata, ngs_annots, igblast_annots, cloneids)
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
        "item",
        "sequence_id_original",
        "partis_clone_id",
        "lineage_group_category",
        "lineage",
        "notes"]
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
    arg("--metadata-isolates", help="CSV with Isolate metadata")
    arg("--metadata-specimens", help="CSV with Specimen metadata")
    arg("--metadata-seqsets", help="CSV with SeqSet metadata")
    arg("-n", "--ngs-annotations", help="optional CSV with Lineage info for known NGS sequences")
    arg("-A", "--igblast-airr", help="optional AIRR tsv.gz from IgBLAST to prefer for annotations")
    arg("-a", "--all", action="store_true",
        help="keep all sequences or only those belonging to clones that also include isolates?")
    args = parser.parse_args()
    partis_seq_lineage_info(
        args.input, args.output,
        args.metadata_isolates, args.metadata_specimens, args.metadata_seqsets,
        args.ngs_annotations, args.igblast_airr, keep_all=args.all)

if __name__ == "__main__":
    main()
