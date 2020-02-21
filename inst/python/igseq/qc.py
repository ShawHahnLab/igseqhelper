"""
Various QC checks.

Currently just some helper functions to summarize the possible behavior of
cutadapt's quality trimming function.

https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm
https://github.com/marcelm/cutadapt/blob/master/src/cutadapt/qualtrim.pyx#L6
"""

import gzip
import csv
import re
from random import random
from Bio import SeqIO
from cutadapt import qualtrim

def make_qualtrim_grid(fqgz_in, qual_breaks=None, len_breaks=None):
    """Tally how many reads would be trimmed to what lengths at what quality cutoffs.

    This produces a grid (in the form of a per-quality-beak dictionary with
    each value being a per-length dictionary containing read counts).  Each
    per-quality-value dictionary contains numbers summing to the total number
    of reads, since all reads are considered for each quality value.
    """

    # Go from highest qual cutoff to lowest
    if not qual_breaks:
        qual_breaks = list(range(0, 42))
    qual_breaks = sorted(qual_breaks)[::-1]
    # Go from lowest to highets length
    if not len_breaks:
        len_breaks = list(range(0, 401))
    len_breaks = sorted(len_breaks)

    tally = {b: {lenb: 0 for lenb in len_breaks} for b in qual_breaks}
    with gzip.open(fqgz_in, "rt") as f_in:
        for record in SeqIO.parse(f_in, "fastq"):
            qual = [chr(val+33) for val in record.letter_annotations["phred_quality"]]
            qual = ''.join(qual)
            for cutoff in qual_breaks:
                _, trim = qualtrim.quality_trim_index(qual, 0, cutoff)
                diffs = [abs(len_break - trim) for len_break in len_breaks]
                len_break = len_breaks[diffs.index(min(diffs))]
                tally[cutoff][len_break] += 1
    return tally

def make_qualtrim_csv(fqgz_in, fp_csv, qual_breaks=None, len_breaks=None):
    """CSV writer for make_qualtrim_grid.

    This writes a CSV file with quality cutoffs on rows, sequence lengths on
    columns, and counts of occurrences in each cell.
    """
    grid = make_qualtrim_grid(fqgz_in, qual_breaks, len_breaks)
    keys_qual = sorted(grid.keys())
    keys_len = sorted(grid[keys_qual[0]].keys())
    with open(fp_csv, "wt") as f_out:
        writer = csv.writer(f_out)
        header = ["Q"] + keys_len
        writer.writerow(header)
        for key in keys_qual:
            row = [key] + [grid[key][key2] for key2 in keys_len]
            writer.writerow(row)

def load_rearr_centroids_by_seq(fp_tsv):
    """Load rearrangements.tsv from SONAR, extract centroid ID and original sequence ID.

    This creates a dictionary of sequence ID -> centroid ID mappings.
    """
    # From a quick check it looks like these two columns put us in the range of
    # 10s of MB of text input, so I think we should be safe to just load it all in.
    # I think.
    with open(fp_tsv) as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        pairs = [(re.sub(r"\|.*$", "", row["source_id"]), row["centroid"]) for row in reader]
    return dict(pairs)

def get_rearr_centroids_by_raw_reads(fp_tsv, fps_fqgz, fp_out):
    """Match raw read IDs for a specimen with SONAR's centroid IDs.

    fp_tsv: SONAR rearrangements.tsv file
    fps_fqgz: list of raw fastq.gz files containing IDs for the original reads
    fp_out: CSV output filename

    This adds entries for read IDs not present in SONAR's set.  We can use this
    to investigate the depth of sampling within and between cells.  This is
    very space-inefficient, but simple to work with.
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
