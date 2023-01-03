"""
Summarizing and reporting helper functions - SONAR.
"""

import csv
import re
import logging
import gzip
import itertools
from statistics import median
from random import random
from Bio import SeqIO
from ..util import parse_seq_desc

LOGGER = logging.getLogger(__name__)

def sonar_island_stats(fp_output, fp_input_island, fp_input_iddiv, fp_input_fasta, extras=None):
    """Condense the full ID/DIV stats to just those for one island and sumamrize across antibodies.

    fp_output: CSV file
    fp_input_island: islandSeqs.txt file for one specimen
    fp_input_iddiv: goodVJ_unique_id-div.tab for one specimen
    fp_input_fasta: islandSeqs.fa file for one specimen
    extras: dictionary of extra constants to add as columns for convenience
            (e.g. specimen, timepoint)

    Creates one CSV file per specimen/lineage/chain combo, taking the median of
    divergence values across antibodies.
    """
    fieldnames = ["sequence_id", "v_gene", "germ_div", "ab_id_min", "ab_id_median", "ab_id_max"]
    if extras:
        fieldnames = list(extras.keys()) + fieldnames
    else:
        extras = {}
    with open(fp_input_island) as f_in:
        ids = [line.strip() for line in f_in]
    # outer dict: seq ID to attributes
    # each inner dict: key/val pairs from sequence descriptions
    with open(fp_input_fasta) as f_in:
        descs = {rec.id: parse_seq_desc(rec.description) for rec in SeqIO.parse(f_in, "fasta")}
    desc_keys = [val.keys() for val in descs.values()]
    desc_keys = sorted(list(set(itertools.chain(*desc_keys))))
    fieldnames += desc_keys
    with open(fp_input_iddiv) as f_in, open(fp_output, "wt") as f_out:
        reader = csv.DictReader(f_in, delimiter="\t")
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in reader:
            if row["sequence_id"] not in ids:
                continue
            keep = lambda key: key not in ["sequence_id", "v_gene", "germ_div"]
            vals = [float(val) for key, val in row.items() if keep(key)]
            if vals:
                ab_min = min(vals)
                ab_med = round(median(vals), 4)
                ab_max = max(vals)
            else:
                ab_min = ''
                ab_med = ''
                ab_max = ''
            row_out = {**extras, **{
                "sequence_id": row["sequence_id"],
                "v_gene": row["v_gene"],
                "germ_div": row["germ_div"],
                "ab_id_min": ab_min,
                "ab_id_median": ab_med,
                "ab_id_max": ab_max}}
            for key in desc_keys:
                row_out[key] = descs[row["sequence_id"]].get(key, "")
            writer.writerow(row_out)

def sonar_island_summary(fp_output_csv, fps_input_csv):
    """Further condense ID/DIV stats to one file per lineage.

    fp_output_csv: CSV for per-lineage ID/DIV information.
    fps_input_csv: list of per-specimen ID/DIV stats.  See sonar_island_stats.

    Each specimen file is collapsed to one row so this summarizies the shift
    across timepoints.
    """
    fieldnames = [
        "specimen", "timepoint", "total",
        "germ_div_min", "germ_div_max", "germ_div_median",
        "ab_id_min", "ab_id_max", "ab_id_median"]
    with open(fp_output_csv, "wt") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        rows_out = []
        for fp_in in fps_input_csv:
            rows_out.append(_sonar_island_summary_row(fp_in))
        rows_out = [row for row in rows_out if row["total"] > 0]
        def sorter(row):
            week = re.search("WK([0-9]+)", row["specimen"])
            if week:
                week = int(week.group(1))
            else:
                week = -1
            return (week, row["specimen"])
        rows_out = sorted(rows_out, key=sorter)
        writer.writerows(rows_out)

def _sonar_island_summary_row(fp_in):
    germ_divs = []
    ab_ids_meds = []
    ab_ids_mins = []
    ab_ids_maxes = []
    specimen = ""
    timepoint = ""
    with open(fp_in) as f_in:
        reader = csv.DictReader(f_in)
        for row in reader:
            germ_divs.append(float(row["germ_div"]))
            ab_ids_meds.append(float(row["ab_id_median"]))
            ab_ids_mins.append(float(row["ab_id_min"]))
            ab_ids_maxes.append(float(row["ab_id_max"]))
            # just take the last specimen and timepoint given (if any) since
            # they should be constant per file
            specimen = row.get("specimen", "")
            timepoint = row.get("timepoint", "")
    if germ_divs:
        row_out = {
            "specimen": specimen,
            "timepoint": timepoint,
            "total": len(germ_divs),
            "germ_div_min": round(min(germ_divs), 4),
            "germ_div_max": round(max(germ_divs), 4),
            "germ_div_median": round(median(germ_divs), 4),
            "ab_id_min": round(min(ab_ids_mins), 4),
            "ab_id_max": round(max(ab_ids_maxes), 4),
            "ab_id_median": round(median(ab_ids_meds), 4)}
    else:
        row_out = {
            "specimen": specimen,
            "timepoint": timepoint,
            "total": 0,
            "germ_div_min": '',
            "germ_div_max": '',
            "germ_div_median": '',
            "ab_id_min": '',
            "ab_id_max": '',
            "ab_id_median": ''}
    return row_out
