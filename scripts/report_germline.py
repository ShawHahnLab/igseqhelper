#!/usr/bin/env python

"""
Tabulate per-subject per-segment germline details.

Make a brief CSV of key/value pairs with information about personalized V+J
germline sequences (e.g. from IgDiscover) versus those in the starting database
those sequences were based on.
"""

import sys
from pathlib import Path
from igseq.record import RecordReader, RecordWriter

def inany(fragment, items):
    """Is the given fragment contained in any of the given items?"""
    return any(fragment in txt for txt in items)

def report_germline(subject, locus, segment, fasta_in, info_csv_in, csv_out):
    """Tabulate per-subject per-segment germline details"""
    fasta_in = Path(fasta_in)
    segment = Path(fasta_in).stem
    seqs = set()
    seqs_trunc = set()
    seqs_orig = set()

    # load all seqs, and gather truncated versions based on segment too
    with RecordReader(fasta_in) as reader:
        for row in reader:
            seqs.add(row["sequence"])
            if segment == "V":
                seqs_trunc.add(row["sequence"][:-10])
            elif segment == "J":
                seqs_trunc.add(row["sequence"][10:])

    # infer reference name and find the starting DB FASTA if possible
    src_path = None
    ref_name = None
    with RecordReader(info_csv_in) as reader:
        for row in reader:
            if row["key"] == segment:
                src_path = Path(row["value"])
                if "oldmethod" in src_path.parts:
                    ref_name = "sonarramesh"
                else:
                    ref_name = src_path.parts[2]
    # if found, gather starting seqs
    if ref_name:
        chain_type = {"IGH": "mu", "IGK": "kappa", "IGL": "lambda"}[locus]
        start_fasta = Path(f"analysis/igdiscover/{ref_name}/{chain_type}/{segment}.fasta")
        if start_fasta.exists():
            with RecordReader(start_fasta) as reader:
                for row in reader:
                    seqs_orig.add(row["sequence"])

    num_seqs_total = len(seqs)
    if seqs_orig:
        num_seqs_ref = len(seqs_orig)
        num_seqs_novel = len(seqs - seqs_orig)
        num_seqs_trunc = len(seqs_trunc)
        fract_seqs_novel = num_seqs_novel/num_seqs_total
        fract_seqs_novel_trunc = sum(
            1 for s in seqs_trunc if not inany(s, seqs_orig))/len(seqs_trunc)
    else:
        num_seqs_ref = None
        num_seqs_novel = None
        fract_seqs_novel = None
        fract_seqs_novel_trunc = None

    # write output
    with RecordWriter(csv_out, "csv") as writer:
        writer.write({"key": "subject", "value": subject})
        writer.write({"key": "locus", "value": locus})
        writer.write({"key": "segment", "value": segment})
        writer.write({"key": "reference", "value": ref_name})
        writer.write({"key": "seqs_total", "value": num_seqs_total})
        writer.write({"key": "seqs_ref", "value": num_seqs_ref})
        writer.write({"key": "seqs_novel", "value": num_seqs_novel})
        writer.write({"key": "seqs_trunc", "value": num_seqs_trunc})
        txt = "" if fract_seqs_novel is None else f"{fract_seqs_novel:0.3f}"
        writer.write({"key": "fract_seqs_novel", "value": txt})
        txt = "" if fract_seqs_novel_trunc is None else f"{fract_seqs_novel_trunc:0.3f}"
        writer.write({"key": "fract_seqs_novel_trunc", "value": txt})

if __name__ == "__main__":
    report_germline(*sys.argv[1:])
