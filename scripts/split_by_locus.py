#!/usr/bin/env python

import argparse
from pathlib import Path
from igseq import igblast
from igseq import vdj
from igseq.record import RecordReader, RecordWriter

def _split_single(input_fqgz, output_handles, vdj_ref_paths, organism, threads):
    with igblast.setup_db_dir(vdj_ref_paths) as (db_dir, _):
        with igblast.run_igblast(
                db_dir, organism, input_fqgz[0], threads, extra_args=["-outfmt", "19"]) as proc:
            with \
                    RecordReader(input_fqgz[0]) as rdr_orig, \
                    RecordReader(proc.stdout, fmt="tsv") as rdr:
                for read, row in zip(rdr_orig, rdr):
                    locus = row["locus"] if row["locus"] in ("IGH", "IGK", "IGL") else "other"
                    output_handles[(locus, "R1")].write(read)
            proc.wait()

def _split_paired(input_fqgz, output_handles, vdj_ref_paths, organism, threads):
    def doigblast(path):
        return igblast.run_igblast(db_dir, organism, path, threads, extra_args=["-outfmt", "19"])
    with igblast.setup_db_dir(vdj_ref_paths) as (db_dir, _):
        with doigblast(input_fqgz[0]) as proc_r1, doigblast(input_fqgz[1]) as proc_r2:
            with \
                    RecordReader(input_fqgz[0]) as rdr_r1_orig, \
                    RecordReader(input_fqgz[1]) as rdr_r2_orig, \
                    RecordReader(proc_r1.stdout, fmt="tsv") as rdr_r1, \
                    RecordReader(proc_r2.stdout, fmt="tsv") as rdr_r2:
                for rd_r1, rd_r2, row_r1, row_r2 in zip(rdr_r1_orig, rdr_r2_orig, rdr_r1, rdr_r2):
                    loci = "/".join([row_r1["locus"], row_r2["locus"]])
                    locus = {
                        "IGH/IGH": "IGH", "/IGH": "IGH",
                        "IGK/IGK": "IGK", "/IGK": "IGK",
                        "IGL/IGL": "IGL", "/IGL": "IGL",
                        }.get(loci, "other")
                    output_handles[(locus, "R1")].write(rd_r1)
                    output_handles[(locus, "R2")].write(rd_r2)
            proc_r1.wait()
            proc_r2.wait()

def _setup_output_paths(input_fqgz, output_pattern=None):
    if len(input_fqgz) == 1:
        rps = ("R1", )
        if not output_pattern:
            output_pattern = "{locus}.fastq.gz"
    elif len(input_fqgz) == 2:
        rps = ("R1", "R2")
        output_pattern = "{locus}.{rp}.fastq.gz"
    else:
        raise ValueError("Input should be single or R1/R2 pair of fastq.gz paths")
    output_paths = {}
    for locus in ["IGH", "IGK", "IGL", "other"]:
        for rp in rps:
            output_paths[(locus, rp)] = Path(output_pattern.format(locus=locus, rp=rp))
    return output_paths

def split_by_locus(input_fqgz, output_pattern=None, species="rhesus", threads=1):
    output_paths = _setup_output_paths(input_fqgz, output_pattern)
    for path in output_paths.values():
        if path.exists():
            raise ValueError(f"Output already exists: {path}")
    ref_paths = [igblast.fuzzy_species_match(species)]
    attrs_list = vdj.parse_vdj_paths(ref_paths)
    vdj_ref_paths = [attrs["path"] for attrs in attrs_list]
    output_handles = {}
    try:
        for key, path in output_paths.items():
            path.parent.mkdir(exist_ok=True, parents=True)
            output_handles[key] = RecordWriter(path)
            output_handles[key].open()
        organism = igblast.detect_organism(None, species)
        splitter = _split_single if len(input_fqgz) == 1 else _split_paired
        splitter(input_fqgz, output_handles, vdj_ref_paths, organism, threads)
    finally:
        for handle in output_handles.values():
            handle.close()

def main():
    """CLI for split_by_locus()"""
    parser = argparse.ArgumentParser()
    addarg = parser.add_argument
    addarg("input", nargs="+", help="Single input fastq.gz or R1 and R2 in order")
    addarg("-o", "--output-pattern",
        help="File path pattern for output (default is "
        "\"{locus}.{rp}.fastq.gz\" for paired and "
        "\"{locus}.fastq.gz\" otherwise")
    addarg("-S", "--species", default="rhesus",
        help="species to use (human or rhesus)")
    addarg("-t", "--threads", type=int, default=1,
        help="number of threads for parallel processing (default: 1)")
    args = parser.parse_args()
    split_by_locus(args.input, args.output_pattern, args.species, args.threads)

if __name__ == "__main__":
    main()
