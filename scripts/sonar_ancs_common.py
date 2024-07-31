#!/usr/bin/env python

"""
Filter inferred ancestors to those leading to a given clade

Assumes SONAR/IgPhyML sequence IDs that contain each subtree's nodes.  The
clade sequences are only used for sequence IDs.
"""

import argparse
from igseq.__main__ import rewrap
from igseq.record import RecordReader, RecordWriter

def __setup_arg_parser():
    parser = argparse.ArgumentParser(
        description=rewrap(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", "--ancestors", help="SONAR/IgPhyML inferred ancestors seqs/rows")
    parser.add_argument("-c", "--clade", help="target clade seqs/rows")
    parser.add_argument("-o", "--output", help="filtered output")
    return parser

def sonar_ancs_common(path_ancestors, path_clade, path_output):
    with RecordReader(path_clade) as reader:
        clade_ids = {rec["sequence_id"] for rec in reader}
    output = []
    with RecordReader(path_ancestors) as reader:
        for rec in reader:
            depth, nodes, dist = rec["sequence_id"].split(";")
            nodes = nodes.split(",")
            if all(clade_id.upper() in nodes for clade_id in clade_ids):
                rec["depth"] = int(depth)
                rec["dist"] = float(dist)
                output.append(rec)
    digits = 0
    if output:
        output = sorted(output, key=lambda rec: rec["depth"])
        digits = len(str(output[-1]["depth"]))
    with RecordWriter(path_output) as writer:
        for rec in output:
            depth = str(rec["depth"]).zfill(digits)
            rec["sequence_id"] = f"ancestor{depth}"
            writer.write(rec)

def main(arglist=None):
    parser = __setup_arg_parser()
    if arglist is None:
        args, args_extra = parser.parse_known_args()
    else:
        args, args_extra = parser.parse_known_args(arglist)
    sonar_ancs_common(args.ancestors, args.clade, args.output)

if __name__ == "__main__":
    main()
