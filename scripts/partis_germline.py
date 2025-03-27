#!/usr/bin/env python

"""Prep partis germline dir from V(D)J FASTA files"""

import argparse
from tempfile import NamedTemporaryFile
from igseq.aux import make_aux_file
from igseq.igblast import setup_db_dir, run_igblast
from igseq.record import RecordReader, RecordWriter
from igseq.vdj import parse_vdj_paths

KEYS = ["gene", "cyst_position", "tryp_position", "phen_position", "aligned_seq"]

def _prep_j(j_fasta):
    # We can use igseq's IgBLAST-related functionality to get the conserved W/F
    # positions in J
    out = []
    with NamedTemporaryFile() as tmp:
        make_aux_file(j_fasta, tmp.name)
        with open(tmp.name, encoding="ASCII") as f_in:
            for line in f_in:
                fields = line.strip().split("\t")
                out_row = {key: "" for key in KEYS}
                # tryptophan for heavy chain, phenylalanine for light chains
                out_aa = {"JH": "tryp", "JK": "phen", "JL": "phen"}[fields[2]]
                out_row.update({"gene": fields[0], f"{out_aa}_position": fields[3]})
                out.append(out_row)
    return out

def _prep_v(threads, v_fasta):
    out = []
    # For the conserved C positions in V, we can IgBLAST and use the AIRR TSV
    # output
    with setup_db_dir([attrs["path"] for attrs in parse_vdj_paths(["rhesus"])]) as (db_dir, _):
        with run_igblast(
                db_dir, "rhesus_monkey", v_fasta,
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

def _copy_segment(src, dir_out, locus, seg):
    locus = locus.lower()
    seg = seg.lower()
    dst = f"{dir_out}/{locus}/{locus}{seg}.fasta"
    # Ignore D if not IGH
    if locus != "igh" and seg == "d":
        return
    with open(src, encoding="ASCII") as f_in, open(dst, "w", encoding="ASCII") as f_out:
        f_out.write(f_in.read())

def partis_germline(locus_inputs, dir_out, threads=1):
    """Prep partis germline dir from V/D/J FASTA files"""
    for locus, seg_paths in locus_inputs.items():
        out = _prep_j(seg_paths["j"])
        out += _prep_v(threads, seg_paths["v"])
        with RecordWriter(f"{dir_out}/{locus}/extras.csv") as writer:
            for row in out:
                writer.write(row)
        for seg, path in seg_paths.items():
            _copy_segment(path, dir_out, locus, seg)

def main():
    """CLI for partis_germline"""
    parser = argparse.ArgumentParser()
    addarg = parser.add_argument
    addarg("locus", type=lambda obj: obj.lower(), choices=["igh", "igk", "igl"],
        help="igh/igk/igl locus")
    addarg("v", help="V FASTA input")
    addarg("d", help="D FASTA input (ignored for light chains)")
    addarg("j", help="J FASTA input")
    addarg("output", help="output directory to create per-locus subdirs under")
    addarg("-t", "--threads", type=int, default=1)
    args = parser.parse_args()
    locus_inputs = {args.locus: {"v": args.v, "d": args.d, "j": args.j}}
    partis_germline(locus_inputs, args.output, args.threads)

if __name__ == "__main__":
    main()
