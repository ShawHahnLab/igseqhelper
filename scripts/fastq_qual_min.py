#!/usr/bin/env python

"""
Keep only records where all quality scores are equal or higher than a given minimum.
"""

import re
import argparse
from igseq.record import RecordReader, RecordWriter
from igseq.util import save_counts, parse_quals

QUALMIN=10

def fastq_qual_min(path_input, path_output, path_counts, qualmin=QUALMIN):
    """Keep only records where all quality scores are equal or higher than a given minimum"""
    cts = {
        "Category": "qualmin",
        "Sample": re.sub(r"\.fastq\.gz", "", re.sub(".*/", "", path_input)),
        "NumSeqs": 0}
    cts_input = {"Item": "input"}
    cts_output = {"Item": "output"}
    cts_input.update(cts)
    cts_output.update(cts)
    with RecordReader(path_input) as reader, RecordWriter(path_output) as writer:
        for record in reader:
            quals = parse_quals(record["sequence_quality"])
            cts_input["NumSeqs"] += 1
            if all(qual >= qualmin for qual in quals):
                cts_output["NumSeqs"] += 1
                writer.write(record)
    save_counts(path_counts, [cts_input, cts_output])

def main():
    """CLI for fastq_qual_min"""
    parser = argparse.ArgumentParser()
    addarg = parser.add_argument
    addarg("input", help="input FASTQ path")
    addarg("output", help="output FASTQ path")
    addarg("-c", "--countsfile", help="file to write read counts to")
    addarg("-Q", "--qualmin", type=int, default=QUALMIN,
        help="minimum quality score required per base for each read")
    args = parser.parse_args()
    fastq_qual_min(args.input, args.output, args.countsfile, args.qualmin)

if __name__ == "__main__":
    main()
