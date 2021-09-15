"""
Summarizing and reporting helper functions - demultiplexing.
"""

import csv
import gzip
from collections import defaultdict

def make_barcode_summary(csv_out, csvs_in, sequences, samples):
    """Combine chunked Seq ID/Fwd Barcode/Rev Barcode tables in one summary.

    csv_out: path to write summary CSV
    csvs_in: list of paths for csv.gz files created by demux.annotate
    sequences: list of dicts of sequence attributes
    samples: list of dicts of sample attributes for this run

    The output CSV has these columns:

      * BCFWD: forward barcode sequence
      * BCREV: reverse barcode sequence
      * BCFWDName: name of forward barcode for recognized barcodes
      * BCREVName: name of reverse barcode for recognized barcodes
      * Total: total read count per barcode pair
      * Sample: sample name for barcode pair for recognized pairs, if samples
                given
    """
    # Tally up the total counts per recognized barcode pair across chunked
    # files
    totals = _load_barcode_totals(csvs_in)
    # Combine into a single per-run summary, including rows for unrecognized
    # cases for fwd/rev/both
    output = _make_barcode_summary(totals, sequences, samples)
    # Write as CSV
    fieldnames = ["BCFWD", "BCREV", "BCFWDName", "BCREVName", "Total", "Sample"]
    with open(csv_out, "wt") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)

def _load_barcode_totals(csvs_in):
    totals = defaultdict(int)
    for csv_in in csvs_in:
        with gzip.open(csv_in, "rt") as f_in:
            reader = csv.DictReader(f_in)
            for row in reader:
                totals[(row["BCFWD"], row["BCREV"])] += 1
    return totals

def _make_barcode_summary(totals, sequences, samples):
    known_fwd = {val["Seq"]: key for key, val in sequences.items() if "BC_" in key}
    known_rev = {val["Seq"]: key for key, val in sequences.items() if "i7_" in key}
    output = []
    def _make_barcode_summary_row(fwd_seq, rev_seq):
        return {
            "BCFWD": fwd_seq,
            "BCREV": rev_seq,
            "BCFWDName": known_fwd.get(fwd_seq),
            "BCREVName": known_rev.get(rev_seq),
            "Total": totals[(fwd_seq, rev_seq)],
            "Sample": _get_samp(known_fwd.get(fwd_seq), known_rev.get(rev_seq), samples)}
    for fwd_seq in known_fwd:
        for rev_seq in known_rev:
            # row for each recognized barocde pair
            output.append(_make_barcode_summary_row(fwd_seq, rev_seq))
        # extra row for each recognized fwd barcode with no recognized rev
        output.append(_make_barcode_summary_row(fwd_seq, ""))
    # extra row for no recognized fwd OR rev barcode
    output.append(_make_barcode_summary_row("", ""))
    # extra rows for no recognized fwd barcode but recognized rev
    for rev_seq in known_rev:
        output.append(_make_barcode_summary_row("", rev_seq))
    return output

def _get_samp(fwd, rev, samples):
    for s_attrs in samples.values():
        if s_attrs["BarcodeFwd"] == fwd and s_attrs["BarcodeRev"] == rev:
            return s_attrs["Sample"]
    return None
