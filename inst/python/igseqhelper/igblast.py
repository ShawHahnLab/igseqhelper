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
