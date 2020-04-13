"""
IgBLAST helpers.
"""

import csv
import re


def split_igblast_seqids_by_chain(input_fp, heavy_out, light_out):
    """Divide Seq IDs from IgBLAST AIRR table into two files."""
    with open(input_fp) as f_in, \
        open(heavy_out, "wt", newline="") as f_out_h, \
        open(light_out, "wt", newline="") as f_out_l:
        reader = csv.DictReader(f_in, delimiter="\t")
        for row in reader:
            seqid = re.sub(r"^reversed\|", "", row["sequence_id"])
            if "IGKV" in row["v_call"] or "IGLV" in row["v_call"]:
                f_out_l.write("%s\n" % seqid)
            elif "IGHV" in row["v_call"]:
                f_out_h.write("%s\n" % seqid)

def tabulate_igblast(input_fp, output_fp):
    """Make CSV from igblastn output with -outfmt 7.

    7 is supposed to be tabular output, but it's more like a giant list of
    tables of varying widths.  Here we'll try to standardize a portion of that
    info into an actual grid.
    """
    fields_out = [
        "Query", "ChainType", "TopV", "TopD", "TopJ",
        "Chain", "StopCodon", "VJFrame", "Productive", "Strand"]
    with open(input_fp) as f_in, open(output_fp, "wt", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fields_out)
        writer.writeheader()
        query = None
        rearr_summary = False
        for line in f_in:
            match = re.match("^# Query: (.*)$", line)
            if match:
                query = match.group(1)
                continue
            match = re.match(
                r"^# V-\(D\)-J rearrangement summary for query sequence \((.*)\).*$", line)
            if match:
                heavy = "Top D" in match.group(1)
                rearr_summary = True
                continue
            if rearr_summary:
                values = line.strip().split("\t")
                fields = [
                    "TopV", "TopD", "TopJ",
                    "Chain", "StopCodon", "VJFrame",
                    "Productive", "Strand"]
                if not heavy:
                    del fields[1]
                row = dict(zip(fields, values))
                if not heavy:
                    row["TopD"] = ""
                row["Query"] = query
                row["ChainType"] = {True: "heavy", False: "light"}[heavy]
                writer.writerow(row)
                rearr_summary = False
