#!/usr/bin/env python

"""
Report sequence and lineage info from partis, merging metadata for our seqs+isolates.

This notes rows that seem like duplicates or conflicting sequences based on
cell barcodes, for 10x, but doesn't actually remove them.
"""

import re
import sys
import gzip
import argparse
from collections import defaultdict
from csv import DictReader, DictWriter

def get_family(txt):
    if not txt:
        return ""
    family = re.match("(IG[HKL][VDJ][0-9]+)", txt)
    return family.group(1) if family else ""

def _load_metadata(csv_path, key=None):
    things = {}
    if csv_path:
        first = lambda row: list(row.keys())[0]
        with open(csv_path, encoding="ASCII") as f_in:
            things = {row[key or first(row)]: row for row in DictReader(f_in)}
    return things

def _load_igblast_airr(airr_path, key="sequence_id"):
    if airr_path:
        with gzip.open(airr_path, "rt", encoding="ASCII") as f_in:
            return {row[key]: row for row in DictReader(f_in, delimiter="\t")}
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
    # will assume any 16 NT motif inside a seq ID, delimited by dashes or
    # underscores, is a 10x cell barcode
    cell_barcode = ""
    if (match := re.search(r"[-_]([ACTG]{16})[-_]", seqid_in)):
        cell_barcode = match.group(1)
    row_out = {
        "sequence_id": seqid_in,
        "sequence_id_original": seqid,
        "cell_barcode": cell_barcode,
        "category": category,
        "item": item,
        "timepoint": attrs.get("Timepoint", timepoint_seqid),
        "notes": [],
        "lineage": "",
        "sequence_light": ""}
    # Catch the special case of seqset entries (10x) that have been added to
    # the isolates table.  Watch out for duplicates though.
    if category == "seqset":
        isolate_alt_name = attrs["Subject"] + "-wk" + attrs["Timepoint"] + "-" + seqid
        isolate_map = {r["AltName"]: r for r in metadata["isolates"].values() if r["AltName"]}
        if isolate_alt_name in isolate_map:
            row_out["notes"].append("Matched to isolate via seqset entry seq ID")
            attrs = isolate_map[isolate_alt_name]
            row_out["category"] = "isolate"
            row_out["item"] = ""
            row_out["sequence_id_original"] = attrs["Isolate"]
            row_out["lineage"] = attrs["Lineage"]
    if attrs.get("Skip") == "TRUE":
        # skip entries if that's noted in their metadata (isolates
        # we don't want to actually analyze, basically); generally
        # they shouldn't get this far anyway, but if so, we'll
        # exclude them now
        return None
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
    # Note paired light seqs where available (this will just be for isolates)
    row_out["sequence_light"] = attrs.get("LightSeq", "")
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
    row_out.update({
        "sequence": igblast_attrs["sequence"],
        "v_family": get_family(igblast_attrs["v_call"]),
        "j_family": get_family(igblast_attrs["j_call"]),
        "v_identity": igblast_attrs["v_identity"],
        "d_call": igblast_attrs["d_call"],
        "junction_aa": igblast_attrs["junction_aa"],
        "junction_aa_length":
            len(igblast_attrs["junction_aa"]) if igblast_attrs["sequence"] else None})

def __include_custom_attrs(row_out, custom_annots):
    # Do we have our own manual annotations for this ID?  If so, take the
    # lineage from there if not otherwise specified
    # (now using the same full, unambiguous ID used in my partis rules, so we
    # can generalize this beyond just the NGS rows)
    if (custom_attrs := custom_annots.get(row_out["sequence_id"])):
        # also sanity check with sequence content if present
        if custom_attrs.get("sequence") and \
                row_out.get("sequence") and \
                custom_attrs["sequence"] not in row_out["sequence"]:
            print(row_out["sequence"])
            print(custom_attrs["sequence"])
            raise ValueError(f"Sequence mismatch for {row_out['sequence_id']}")
        if not row_out["lineage"]:
            row_out["lineage"] = custom_attrs.get("Lineage", "")

def __include_isolate_light_attrs(row_out, isolate_light_annots):
    attrs = isolate_light_annots.get(row_out["sequence_light"], {})
    row_out.update({
    "light_v_family": get_family(attrs.get("v_call")),
    "light_j_family": get_family(attrs.get("j_call")),
    "light_v_identity": attrs.get("v_identity"),
    "light_junction_aa": attrs.get("junction_aa"),
    "light_junction_aa_length":
        len(attrs.get("junction_aa", "")) if row_out["sequence_light"] else None})

def _check_for_duplicated_isolates(out):
    # Sanity-check the isolates to ensure we don't have duplicates.  I worry
    # this could happen for the 10x sequences that could be present both in the
    # "seqsets" files and also stored as isolates.
    isolate_tally = defaultdict(int)
    for row in out:
        if row["category"] == "isolate":
            isolate_tally[row["sequence_id_original"]] += 1
    isolate_tally = {key: val for key, val in isolate_tally.items() if val > 1}
    if isolate_tally:
        print("Duplicated isolates in output!")
        for isolate, num in isolate_tally.items():
            print(f"  {isolate}: {num}")

def _exclude_based_on_cell_barcodes(out):
    # for seqset rows with cell barcodes inferred, exclude duplicates, but
    # exclude all rows for the cell if the heavy chain sequences clash.
    exclude_extras = set()
    exclude_clashes = set()
    # first, group those with barcodes, by barcodes
    by_barcode = defaultdict(list)
    for row in out:
        if row["cell_barcode"] and row["category"] == "seqset_10x":
            by_barcode[row["cell_barcode"]].append(row)
    # Confirm only the item identifier and associated long seq ID differ, and
    # if so, mark all but the first for removal (...and excluding notes since I
    # make that a list object)
    skips = ("sequence_id", "item", "notes")
    for chunk in by_barcode.values():
        check = {tuple(((k, v) for k, v in row.items() if k not in skips)) for row in chunk}
        if len(check) != 1:
            # If there's more than one unique case, after excluding the keys
            # that actually should differ, categorize this as a clash between
            # distinct heavy chains for all rows.
            for row in chunk:
                exclude_clashes.add(row["sequence_id"])
        else:
            # Otherwise, just note the extra rows past the first one as
            # duplicates.
            for row in chunk[1:]:
                exclude_extras.add(row["sequence_id"])
    # Mark those cases for exclusion.  (Looping over all rows but this will
    # only set a non-empty string for the applicable 10x cases.)
    if exclude_clashes:
        sys.stderr.write(
            f"Excluding {len(exclude_clashes)} sequences "
            "with mismatched heavy chains within cells\n")
    for row in out:
        row["exclusion_reason"] = ""
        if row["sequence_id"] in exclude_extras:
            row["exclusion_reason"] = "duplicate"
        if row["sequence_id"] in exclude_clashes:
            row["exclusion_reason"] = "clash"

def _prep_seq_lineage_info(clones, metadata, custom_annots, igblast, isolate_light_annots, cloneids):
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
                # skip this one entirely, based on the supplied metadata.
                row_out = __infer_basics_from_metadata(row["sequence_id"], metadata)
                if row_out is None:
                    continue
                # Add additional information with the help of IgBLAST output
                # and (if applicable) manually-defined info on sequences
                __include_igblast_attrs(row_out, igblast)
                __include_custom_attrs(row_out, custom_annots)
                __include_isolate_light_attrs(row_out, isolate_light_annots)
                row_out.update({
                    "partis_clone_id": row["clone_id"] or "",
                    "lineage": row_out["lineage"] or ""})
                out.append(row_out)
    return out

def _assign_lineage_groups(out, lin_prefix=None, auto_group_for=None):
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
                if lin_prefix:
                    row["lineage_group"] = lin_prefix + "-" + row["lineage_group"]
                row["lineage_group_category"] = "automatic"
        elif "" in lins or (auto_group_for and row["category"] in auto_group_for):
            # (Bypassing the empty-lineage-assignment check for given category
            # (originally 10x isolates specifically), so I can always check the
            # partis groupings in case partis merged some of those together.
            # That's only because I've noted a bunch of lineage labels from
            # elsewhere that I haven't yet confirmed.)
            # Otherwise, if it's a mix of assigned and blank lineages, use this
            # clone ID to group by all observed lineage names for this clone.
            # (TODO Wait, should I also span across other clone IDs that
            # overlap by lineage name, too?  Yeah probably.  Currently this is
            # set up so that we could end with some things labeled "linA" and
            # others "linA/linB".  But good enough for now.)
            lineages = set(clone_lineages[row["partis_clone_id"]]) - {""}
            lineages = sorted(lineages)
            # edge case for a few Duke entries that would otherwise be like
            # "DI57-Duke-H035106-K028816/DI57-Duke-H035106-L028115"
            # will instead be like
            # "DI57-Duke-035106"
            duke_pattern = r"([A-Z0-9]+-Duke-H?[0-9]+)-[KL][0-9]*$"
            duke_pattern2 = r"([A-Z0-9]+-Duke-clone[0-9]+)-[KL]$"
            row["lineage_group"] = "/".join(lineages)
            if len(lineages) > 1:
                print(lineages)
                if all(re.match(duke_pattern, lin) for lin in lineages):
                    prefix = re.match(duke_pattern, lineages[0]).group(1)
                    if all(lin.startswith(prefix) for lin in lineages):
                        # if all start like that, then use short form
                        row["lineage_group"] = re.sub(r"-H?([0-9]+)$", r"-\1", prefix)
                elif all(re.match(duke_pattern2, lin) for lin in lineages):
                    # or this
                    # "DH17-Duke-clone103-K/DH17-Duke-clone103-L"
                    # will instead be like
                    # "DH17-Duke-clone103"
                    prefix = re.match(duke_pattern2, lineages[0]).group(1)
                    if all(lin.startswith(prefix) for lin in lineages):
                        # if all start like that, then use short form
                        row["lineage_group"] = prefix
            row["lineage_group_category"] = "partis-grouped"
        else:
            # otherwise just use this row's one lineage as its group name,
            # ignoring partis' grouping
            row["lineage_group"] = row["lineage"]
            row["lineage_group_category"] = "manual"

def _note_uca_diffs(out, uca_annots):
    for row in out:
        row["uca_sequence_id"] = ""
        row["uca_sequence"] = ""
        row["uca_junction_aa_length_diff"] = None
        # Prefer "UCA", then "RUA"
        keys = [row["lineage_group"] + f"_{suf}" for suf in ("UCA", "UCA_Draft", "RUA")]
        for key in keys:
            if (attrs := uca_annots.get(key)):
                len_uca = int(attrs["junction_aa_length"])
                len_ab = int(row["junction_aa_length"])
                diff = len_ab - len_uca
                diff = f"{diff:+}"
                row["uca_sequence_id"] = attrs["sequence_id"]
                row["uca_sequence"] = attrs["sequence"].replace("-", "")
                row["uca_junction_aa_length_diff"] = diff
                break

def _exclude_duke_pair_edge_cases(out):
    # For instances where I have multiple antibody entries with the exact same
    # heavy chain noted (not the same sequence content but literally the same
    # observed heavy chain represented in more than one ab entry), as happens
    # when Duke reported things like IGH+IGK+IGL, ensure I don't include the
    # same one multiple times for a single lineage group.
    # These are isolates with names like "CE89-Duke-H033694-K028155".
    # If this issue comes up, and one of the duplicated cases has a light chain
    # matching the lineage for other members, prefer that one.
    def getlightlocus(row):
        try:
            return row["light_v_family"][2] # "K" or "L"
        except IndexError:
            return None
    chunks = defaultdict(list)
    for row in out:
        chunks[row["lineage_group"]].append(row)
    for rows in chunks.values():
        duke_isol = defaultdict(list)
        lights = defaultdict(int)
        for row in rows:
            if row["category"] == "isolate" and \
                    (match := re.match(r".*-Duke-H([0-9]+)-[KL][0-9]+", row["sequence_id"])):
                # for the Duke ones, group them by heavy ID
                heavy_num = match.group(1)
                duke_isol[heavy_num].append(row)
            else:
                # for the others, just note what light loci are present, if any
                if (light_locus := getlightlocus(row)):
                    lights[light_locus] += 1
        for rows in duke_isol.values():
            # for cases with more than one row for a Duke ab number, keep at
            # most one row
            if len(rows) > 1:
                # is there a clear winner among light loci for the whole lineage?
                # ("K" or "L" or "")
                # if so, prefer that entry
                light_here = sorted(((val, key) for val, key in lights.items()), reverse=True)
                light_here = [pair for pair in light_here if pair[0] == light_here[0][0]]
                light_here = light_here[0][1] if len(light_here) == 1 else None
                rows = sorted(rows, key=lambda row, x=light_here: getlightlocus(row) != x)
                # In any case just take the first one, once we're done with any sorting
                seqid = rows[0]["sequence_id"]
                if not light_here:
                    rows[0]["notes"].append("arbitrarily selected "
                        f"IG{getlightlocus(rows[0])} ab among duplicates")
                for row in rows[1:]:
                    row["exclusion_reason"] = ("duplicate heavy entry with "
                        f"alternate light chain compared with {seqid}")

def _finalize(out):
    for row in out:
        row["notes"] = "; ".join(row["notes"])
    out.sort(key = lambda row: (
        row["lineage_group"] == "",
        row["partis_clone_id"] == "",
        row["lineage_group"],
        row["lineage"],
        row["partis_clone_id"],
        row["sequence_id"]))

def partis_seq_lineage_info(
        airr_in, csv_out,
        metadata_isolates=None, metadata_specimens=None, metadata_seqsets=None,
        csv_custom_annots=None,
        airr_in_igblast=None, airr_in_isolate_light=None, airr_in_uca=None,
        *, lin_prefix=None, keep_all=False, auto_group_for=None):
    """Report sequences with partis clones overlapping with our isolates"""
    # name -> attrs
    metadata = {
        "isolates": _load_metadata(metadata_isolates),
        "specimens": _load_metadata(metadata_specimens),
        "seqsets": _load_metadata(metadata_seqsets),
        }
    # seq ID here -> custom attrs incl. Lineage
    custom_annots = _load_metadata(csv_custom_annots, "sequence_id")
    # seq ID here -> IgBLAST attrs
    igblast_annots = _load_igblast_airr(airr_in_igblast)
    # unique isolate light seq -> AIRR attrs
    # (This is purely to get some very basic attributes for the light chains so
    # we'll just track via sequence)
    isolate_light_annots = _load_igblast_airr(airr_in_isolate_light, "sequence")
    # UCA/RUA seq ID -> AIRR attrs
    # if provided, will note differences (just CDRH3 AA len currently) to UCA/RUA
    uca_annots = _load_igblast_airr(airr_in_uca)
    clones, cloneids = _load_clones_from_partis_airr(airr_in, metadata, keep_all)
    out = _prep_seq_lineage_info(
            clones, metadata, custom_annots, igblast_annots, isolate_light_annots, cloneids)
    _exclude_based_on_cell_barcodes(out)
    _check_for_duplicated_isolates(out)
    _assign_lineage_groups(out, lin_prefix, auto_group_for)
    _note_uca_diffs(out, uca_annots)
    _exclude_duke_pair_edge_cases(out)
    _finalize(out)
    keys_by_chain = ["v_family", "j_family", "v_identity", "junction_aa", "junction_aa_length"]
    keys = [
        "lineage_group",
        "sequence_id",
        "sequence",
        "sequence_light"] + \
            keys_by_chain + \
            [f"light_{key}" for key in keys_by_chain] + [
        "d_call",
        "timepoint",
        "category",
        "item",
        "sequence_id_original",
        "cell_barcode",
        "partis_clone_id",
        "lineage_group_category",
        "lineage",
        "uca_sequence_id",
        "uca_sequence",
        "uca_junction_aa_length_diff",
        "exclusion_reason",
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
    arg("-n", "--custom-annotations", help="optional CSV with Lineage info for known sequences")
    arg("-A", "--igblast-airr", help="optional AIRR tsv.gz from IgBLAST to prefer for annotations")
    arg("-L", "--isolate-light-airr", help="optional AIRR tsv.gz for isolate light chain sequences")
    arg("-U", "--uca-heavy-airr", help="optional AIRR tsv.gz for UCA heavy chain sequences")
    arg("-P", "--lineage-name-prefix",
        help="optional prefix for auto-named lineage groups (like subject ID)")
    arg("-X", "--auto-group-for", nargs="+",
        help="category label(s) to allow merging lineage groups" \
        " even if all have assigned lineages already (e.g. isolate_10x)")
    arg("-a", "--all", action="store_true",
        help="keep all sequences or only those belonging to clones that also include isolates?")
    args = parser.parse_args()
    partis_seq_lineage_info(
        args.input, args.output,
        args.metadata_isolates, args.metadata_specimens, args.metadata_seqsets,
        args.custom_annotations, args.igblast_airr, args.isolate_light_airr, args.uca_heavy_airr,
        lin_prefix=args.lineage_name_prefix, keep_all=args.all, auto_group_for=args.auto_group_for)

if __name__ == "__main__":
    main()
