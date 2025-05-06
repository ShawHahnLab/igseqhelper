#!/usr/bin/env python

"""
Create one .fastq.gz file with raw reads for each SONAR cluster in a rearrangements TSV.

This supports custom processing on the raw data upstream of each of SONAR's
clusters of similar sequences.
"""

import sys
import csv
from pathlib import Path
from collections import defaultdict
from igseq.record import RecordReader, RecordWriter

def _load_clusters(airr_in):
    clusters = defaultdict(list)
    with open(airr_in, encoding="UTF8") as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        for row in reader:
            clusters[row["centroid"]].append(row)
    return clusters

def _load_fqgz_recs(fqgz_paths, clusters):
    # SONAR only keeps one seq ID for each set of identical reads, so we'll
    # structure this as a list of matching reads by sequence content for each
    # cluster.
    fqgz_recs = defaultdict(list)
    for path in fqgz_paths:
        keeps = {}
        for rows in clusters.values():
            for row in rows:
                keeps[row["sequence"]] = row["source_id"]
        with RecordReader(path) as reader:
            for row in reader:
                if row["sequence"] in keeps:
                    readid = keeps[row["sequence"]]
                    fqgz_recs[readid].append(row)
    return fqgz_recs

def sonar_cluster_fastqs(airr_in, dir_out, fqgz_dir_in):
    """Create one .fastq.gz file with raw reads for each SONAR cluster in a rearrangements TSV."""
    dir_out = Path(dir_out)
    clusters = _load_clusters(airr_in)
    fqgzs = Path(fqgz_dir_in).glob("*.fastq.gz")
    fqgz_recs = _load_fqgz_recs(fqgzs, clusters)
    dir_out.mkdir(parents=True, exist_ok=True)
    for centroid, rows in clusters.items():
        qual_recs = []
        for row in rows:
            qual_recs_here = fqgz_recs[row["source_id"]]
            for rec in qual_recs_here:
                # note the cluster and the representative read ID in the
                # description
                rec["sequence_description"] = row["sequence_id"] + " " + row["source_id"]
            qual_recs += qual_recs_here
        with RecordWriter(dir_out / f"{centroid}.fastq.gz") as writer:
            for qual_rec in qual_recs:
                writer.write(qual_rec)

def main():
    """CLI for sonar_cluster_fastqs"""
    sonar_cluster_fastqs(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == "__main__":
    main()
