"""
Summarizing and reporting helper functions - SONAR.
"""

import csv
import re
import logging
import gzip
from random import random
from Bio import SeqIO

LOGGER = logging.getLogger(__name__)

def load_rearr_centroids_by_seq(fp_tsv):
    """Load rearrangements.tsv from SONAR, extract centroid ID and original sequence ID.

    This creates a dictionary of sequence ID -> centroid ID mappings.  Only
    centroids that have an integer cluster_count listed will be included, so
    that we leave out one-offs that have no duplicates or clustered relatives.
    """
    with open(fp_tsv) as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        centroids_keep = {row["centroid"] for row in reader if row["cluster_count"]}
    # From a quick check it looks like these two columns put us in the range of
    # 10s of MB of text input, so I think we should be safe to just load it all in.
    # I think.
    # (This regular expression will strip off the extra stuff in the seq ID from pRESTO.)
    fmt = lambda row: (re.sub(r"\|.*$", "", row["source_id"]), row["centroid"])
    with open(fp_tsv) as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        pairs = [fmt(row) for row in reader if row["centroid"] in centroids_keep]
    return dict(pairs)

def get_rearr_centroids_by_raw_reads(fp_tsv, fps_fqgz, fp_out):
    """Match raw read IDs for a specimen with SONAR's centroid IDs, for centroids with clusters.

    fp_tsv: SONAR rearrangements.tsv file
    fps_fqgz: list of raw fastq.gz files containing IDs for the original reads
    fp_out: CSV output filename

    This adds entries for read IDs not present in SONAR's set.  We can use this
    to investigate the depth of sampling within and between cells.  This is
    very space-inefficient, but simple to work with.  Only centroids that have
    an integer cluster_count listed will be included, so that we leave out
    one-offs that have no duplicates or clustered relatives.
    """
    # this gets most of the way there, but only for sequences SONAR knows about.
    centroids_by_seq = load_rearr_centroids_by_seq(fp_tsv)
    with open(fp_out, "wt") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=["ReadID", "CentroidID"])
        writer.writeheader()
        for fqgz in fps_fqgz:
            with gzip.open(fqgz, "rt") as f_in:
                for record in SeqIO.parse(f_in, "fastq"):
                    centroid = centroids_by_seq.get(record.id, "")
                    writer.writerow({"ReadID": record.id, "CentroidID": centroid})

def rarefy_sonar_clusters(fp_in, fp_out, num_iters=3, step=100000):
    """Subsample ReadID/CentroidID pairs and count centroids observed.

    fp_in: CSV from get_rearr_centroids_by_raw_reads
    fp_out: CSV to write to with results summary
    num_iters: how many times should we randomly sample to a given depth?

    This builds up a rarefaction curve for how many unique SONAR sequence
    centroids we collect for a given number of reads.
    """
    # A simple line count to get the number of reads total
    reads_total = 0
    with open(fp_in) as f_in:
        for _ in f_in:
            reads_total += 1
    # Each subsampling will give us one value for the number of unique
    # centroids encountered.
    with open(fp_out, "wt", buffering=1) as f_out:
        writer = csv.DictWriter(
            f_out,
            fieldnames=["Depth", "ReadsUsed", "Iteration", "UniqueCentroids"])
        writer.writeheader()
        read_num_breaks = list(range(step, reads_total, step)) + [reads_total]
        _rarefy_entry(read_num_breaks, reads_total, num_iters, fp_in, writer)

def _rarefy_entry(read_num_breaks, reads, num_iters, fp_in, writer):
    for read_num in read_num_breaks:
        for attempt in range(num_iters):
            centroids = set()
            read_actual = 0
            with open(fp_in) as f_in:
                reader = csv.reader(f_in)
                for row in reader:
                    # this will be a little fuzzy since we might not sample
                    # exactly read_num, but it should be close enough.
                    if random() < read_num*1.0/reads:
                        read_actual += 1
                        if row[1]:
                            centroids.add(row[1])
            writer.writerow({
                "Depth": read_num,
                "ReadsUsed": read_actual,
                "Iteration": attempt,
                "UniqueCentroids": len(centroids)})
