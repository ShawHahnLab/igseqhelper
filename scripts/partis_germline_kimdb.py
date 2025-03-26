#!/usr/bin/env python

"""Prep partis germline dir from KIMDB sequences"""

import argparse
from igseq.aux import make_aux_file
from igseq.util import DATA
from igseq.igblast import setup_db_dir, run_igblast
from igseq.record import RecordReader, RecordWriter
from igseq.vdj import parse_vdj_paths

KEYS = ["gene", "cyst_position", "tryp_position", "phen_position", "aligned_seq"]

def _prep_j(dir_out):
    # We can use igseq's IgBLAST-related functionality to get the W
    # positions in J
    out = []
    make_aux_file(
        DATA/"germ/rhesus/kimdb/IGH/IGHJ.fasta",
        f"{dir_out}/gl.aux")
    with open(f"{dir_out}/gl.aux", encoding="ASCII") as f_in:
        for line in f_in:
            fields = line.strip().split("\t")
            out_row = {key: "" for key in KEYS}
            out_row.update({"gene": fields[0], "tryp_position": fields[3]})
            out.append(out_row)
    return out

def _prep_v(threads):
    out = []
    # And for the C positions in V, we can IgBLAST and use the AIRR TSV output
    with setup_db_dir([attrs["path"] for attrs in parse_vdj_paths(["kimdb"])]) as (db_dir, _):
        with run_igblast(
                db_dir, "rhesus_monkey", DATA/"germ/rhesus/kimdb/IGH/IGHV.fasta",
                threads, extra_args=["-outfmt", "19"]) as proc:
            with RecordReader(proc.stdout, fmt="tsv") as rdr:
                for row in rdr:
                    # partis denotes cyst_position as the first NT of the
                    # codon that codes for the W, 0-indexed, while fwr3_end
                    # here is the last NT of that codon, 1-indexed.  All
                    # together that works out to 3 less than fwr3_end.
                    pos = int(row["fwr3_end"]) - 3
                    out_row = {key: "" for key in KEYS}
                    out_row.update({"gene": row["sequence_id"], "cyst_position": pos})
                    out.append(out_row)
            proc.wait()
    return out

def partis_germline_kimdb(dir_out, threads=1):
    """Prep partis germline dir from KIMDB sequences"""
    out = _prep_j(dir_out)
    out += _prep_v(threads)
    with RecordWriter(f"{dir_out}/extras.csv") as writer:
        for row in out:
            writer.write(row)
    for key in ("v", "d", "j"):
        src = DATA/f"germ/rhesus/kimdb/IGH/IGH{key.upper()}.fasta"
        dst = f"{dir_out}/igh{key}.fasta"
        with open(src, encoding="ASCII") as f_in, open(dst, "w", encoding="ASCII") as f_out:
            f_out.write(f_in.read())

def main():
    """CLI for partis_germline_kimdb"""
    parser = argparse.ArgumentParser()
    addarg = parser.add_argument
    addarg("output")
    addarg("-t", "--threads", type=int, default=1)
    args = parser.parse_args()
    partis_germline_kimdb(args.output, args.threads)

if __name__ == "__main__":
    main()
