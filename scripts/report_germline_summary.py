#!/usr/bin/env python

"""
Condense per-subject per-segment germline details to one CSV table

Each row is one subject.  Columns are named according to the subject parsed
from each input file, the locus (H, K, L) and segment (V, J):

    {subject}
    {locus}r         Name of reference used for starting DB
    {locus}{seg}rs   Sequences in reference DB
    {locus}{seg}s    Total sequences in personalized germline
    {locus}{seg}ns   ...number that are not in ref DB
    {locus}{seg}fns  ...fraction of total not in ref DB
"""

import sys
from pathlib import Path
from collections import defaultdict
from igseq.record import RecordReader, RecordWriter

LOCI = ("H", "K", "L")
SEGMENTS = ("V", "J")

def _load_info_from_csvs(dirpath):
    info = defaultdict(list)
    for path in Path(dirpath).iterdir():
        if path.is_dir():
            for path2 in path.glob("*.csv"):
                with RecordReader(path2) as reader:
                    row = {row["key"]: row["value"] for row in reader}
                    info[row["subject"]].append(row)
    return info

def _group_by_segment(rows):
    by_segment = {}
    for locus in LOCI:
        for segment in SEGMENTS:
            key = f"{locus}{segment}"
            for row in rows:
                if row["locus"] == f"IG{locus}" and row["segment"] == segment:
                    if key in by_segment:
                        raise ValueError
                    by_segment[key] = row
    return by_segment

def _flatten_by_subject(info):
    out = []
    for subject, rows in info.items():
        row_out = {"subject": subject}
        by_segment = _group_by_segment(rows)
        for locus in LOCI:
            ref = {r["reference"] for k, r in by_segment.items() if k.startswith(f"{locus}")}
            if len(ref) > 1:
                ref = "/".join(ref) + " (???)"
            elif len(ref) == 1:
                ref = ref.pop()
            else:
                ref = None
            row_out[f"{locus}ref"] = ref
            for segment in SEGMENTS:
                key = f"{locus}{segment}"
                row = by_segment.get(key, {})
                row_out.update({
                   f"{key}rs": row.get("seqs_ref"),
                   f"{key}s": row.get("seqs_total"),
                   f"{key}ns": row.get("seqs_novel"),
                   f"{key}fns": row.get("fract_seqs_novel")})
        out.append(row_out)
    return out

def _sort(out):
    for row in out:
        row["_loci"] = sum(bool(row[f"{locus}Vs"]) for locus in LOCI)
        row["_refs"] = "/".join([row[f"{locus}ref"] or "" for locus in LOCI])
        row["_vseqs"] = sum(int(row[f"{locus}Vs"] or 0) for locus in LOCI)
    out.sort(key=lambda row: [row[k] for k in row if k.startswith("_")])
    for row in out:
        keys_rm = [key for key in row if key.startswith("_")]
        for key in keys_rm:
            del row[key]

def report_germline_summary(dirpath, csv_out):
    """Condense per-subject per-segment germline details to one per-subject table"""
    info = _load_info_from_csvs(dirpath)
    out = _flatten_by_subject(info)
    _sort(out)
    with RecordWriter(csv_out, "csv") as writer:
        for row in out:
            writer.write(row)

def main():
    """CLI for report_germline_summary"""
    report_germline_summary(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
