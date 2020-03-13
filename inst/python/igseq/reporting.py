"""
Summarizing and reporting helper functions.
"""

import csv
import re
import logging
import gzip
from random import random
from tempfile import NamedTemporaryFile
from snakemake.shell import shell
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from cutadapt import qualtrim
from igseq.data import load_csv, get_samples_per_run, MetadataError

LOGGER = logging.getLogger(__name__)

def make_qualtrim_grid(fqgz_in, qual_breaks=None, len_breaks=None):
    """Tally how many reads would be trimmed to what lengths at what quality cutoffs.

    This produces a grid (in the form of a per-quality-beak dictionary with
    each value being a per-length dictionary containing read counts).  Each
    per-quality-value dictionary contains numbers summing to the total number
    of reads, since all reads are considered for each quality value.

    https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm
    https://github.com/marcelm/cutadapt/blob/master/src/cutadapt/qualtrim.pyx#L6
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

def counts_sample_summary(
        file_counts_in, counts_sample_summary_out, samples):
    """Take a list of sequence counts for fastq.gz files and summarize per-sample."""
    samples_per_run = get_samples_per_run(samples)
    counts = load_csv(file_counts_in)
    rows = []
    for runid in samples_per_run:
        for sample in samples_per_run[runid] + ["unassigned"]:
            key = "counts/demux/{run}/{sample}.I1.fastq.gz.counts".format(
                run=runid,
                sample=sample)
            cts = int(counts[key]["NumSequences"])
            try:
                cellcount = samples[sample]["SpecimenAttrs"]["CellCount"]
                ratio = divide(cts, cellcount)
            except KeyError:
                cellcount = ""
                ratio = ""
            rows.append({
                "Run": runid,
                "Sample": sample,
                "CellCount": cellcount,
                "NumSequences": cts,
                "Ratio": ratio})
    with open(counts_sample_summary_out, "wt") as f_out:
        writer = csv.DictWriter(
            f_out,
            fieldnames=["Run", "Sample", "NumSequences", "CellCount", "Ratio"])
        writer.writeheader()
        writer.writerows(rows)

def counts_run_summary(counts_sample_summary_in, counts_run_summary_out):
    """Take a per-sample summary of sequence counts and make a per-run version.

    counts_sample_summary_in: CSV file path, as from counts_sample_summary
    counts_run_sumary_out: CSV file path for output
    """
    with open(counts_sample_summary_in) as f_in:
        reader = csv.DictReader(f_in)
        counts_by_sample = list(reader)
    counts_by_run = _get_counts_by_run(counts_by_sample)
    rows = _tally_counts_by_run(counts_by_run)
    with open(counts_run_summary_out, "wt") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["Run", "UnassignedSeqs", "SampleSeqs", "TotalSeqs", "Ratio"])
        writer.writerows(rows)

def _get_counts_by_run(counts_by_sample):
    counts_by_run = {}
    for row_in in counts_by_sample:
        runid = row_in["Run"]
        try:
            cts = counts_by_run[runid]
        except KeyError:
            counts_by_run[runid] = {"samples": 0, "unassigned": 0}
            cts = counts_by_run[runid]
        if row_in["Sample"] == "unassigned":
            key = "unassigned"
        else:
            key = "samples"
        cts[key] += int(row_in["NumSequences"])
    return counts_by_run

def _tally_counts_by_run(counts_by_run):
    rows = []
    for run_name, run_attrs in counts_by_run.items():
        try:
            ratio = divide(run_attrs["unassigned"], run_attrs["samples"])
        except KeyError:
            ratio = ""
        rows.append([
            run_name,
            run_attrs.get("unassigned", ""),
            run_attrs.get("samples", ""),
            run_attrs.get("unassigned", "") + run_attrs.get("samples", ""),
            ratio])
    return rows

def amplicon_summary(input_fps, specimens, regex):
    """Take a list of per-specimen read count files and make a list of summary dictionaries."""
    LOGGER.debug("amplicon_summary: # input_fps: %d", len(input_fps))
    LOGGER.debug("amplicon_summary: # specimens: %d", len(specimens))
    LOGGER.debug("amplicon_summary: regex: %s", regex)
    counts = []
    for fp_in in input_fps:
        match = re.match(regex, fp_in)
        if match:
            chain = match.group(1)
            chain_type = match.group(2)
            specimen = match.group(3)
        else:
            raise ValueError("Couldn't parse specimen details from %s" % fp_in)
        if not specimen in specimens:
            raise MetadataError("Unknown specimen %s" % specimen)
        with open(fp_in) as f_in:
            cts = int(next(f_in))
        counts.append({
            "Subject": specimens[specimen]["Subject"],
            "Timepoint": int(specimens[specimen]["Timepoint"]),
            "Chain": chain,
            "ChainType": chain_type,
            "Specimen": specimen,
            "Seqs": cts})
    fieldnames = [
        "Subject", "Timepoint", "Chain", "ChainType", "Specimen", "Seqs"]
    counts = sorted(counts, key=lambda x: [x[key] for key in fieldnames])
    return counts, fieldnames

def counts_specimen_summary(input_fps, output_fp, specimens):
    """Take a list of per-specimen-R1 read count files and make a summary table."""
    counts, fieldnames = amplicon_summary(
        input_fps,
        specimens,
        r"counts/presto/data/([a-z]+)\.([a-z]+)/([A-Za-z0-9]+).R1.fastq.counts")
    # Add some additional columns of interest for this specific case
    for cts in counts:
        cts["CellCount"] = specimens[cts["Specimen"]]["CellCount"]
        cts["Ratio"] = divide(cts["Seqs"], specimens[cts["Specimen"]]["CellCount"])
    fieldnames.extend(["CellCount", "Ratio"])
    with open(output_fp, "wt") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(counts)

def counts_assembly_summary(input_fps, output_fp, specimens):
    """Take a list of per-specimen assembled sequence count files and make a summary table."""
    counts, fieldnames = amplicon_summary(
        input_fps,
        specimens,
        r"counts/presto/assemble/([a-z]+)\.([a-z]+)/([A-Za-z0-9]+)_assemble-pass.fastq.counts")
    with open(output_fp, "wt") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(counts)

def counts_presto_qual_summary(input_fps, output_fp, specimens):
    """Take a list of per-specimen presto quality sequence count files and make a summary table."""
    counts, fieldnames = amplicon_summary(
        input_fps,
        specimens,
        r"counts/presto/qual/([a-z]+)\.([a-z]+)/([A-Za-z0-9]+)_quality-pass.fastq.counts")
    with open(output_fp, "wt") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(counts)

def divide(val1, val2, fmt="{:.6f}"):
    """Divide val1 by val2 as floats and return formatted string."""
    num = float(val1)/float(val2)
    num = fmt.format(num)
    return num

def gather_antibodies(
        subject, chain, antibody_isolates, output_fasta):
    """Gather antibody sequences from metadata for use in alignments."""
    #lineages = {k: val for k, val in antibody_lineages.items() if val["Subject"] == subject}
    #in_lineages = lambda v: v["AntibodyLineage"] in lineages.keys()
    #isolates = {k: val for k, val in antibody_isolates.items() if in_lineages(val)}
    is_subject = lambda val: val["AntibodyLineageAttrs"]["Subject"] == subject
    isolates = {k: val for k, val in antibody_isolates.items() if is_subject(val)}
    lineages = {val["AntibodyLineage"]: val["AntibodyLineageAttrs"] for val in isolates.values()}
    if chain == "heavy":
        seq_col = "HeavySeq"
        cons_col = "HeavyConsensus"
    else:
        seq_col = "LightSeq"
        cons_col = "LightConsensus"
    with open(output_fasta, "wt") as f_out:
        for lineage_name, lineage_attrs in lineages.items():
            record = SeqRecord(Seq(lineage_attrs[cons_col]), id=lineage_name, description="")
            SeqIO.write(record, f_out, "fasta")
        for isolate_name, isolate_attrs in isolates.items():
            record = SeqRecord(Seq(isolate_attrs[seq_col]), id=isolate_name, description="")
            SeqIO.write(record, f_out, "fasta")

def align_next_segment(antibodies_fasta, aligned_fasta, alleles_fasta, output_fasta):
    """Align discovered alleles from a segment to the regions not yet aligned to.

    antibodies_fasta: path to FASTA for known antibody sequences
    aligned_fasta: path to FASTA for previous segment's alignment
    alleles_fasta: path to FASTA for discovered alleleles for next segment
    output_fasta: path to FASTA to save new alignment to
    """
    antibodies = list(SeqIO.parse(antibodies_fasta, "fasta"))
    alignment = list(SeqIO.parse(aligned_fasta, "fasta"))
    antibody_ids = {record.id for record in antibodies}

    # Find the farthest-left right edge of the previously-aligned set of alleles.
    right_ends = []
    for record in alignment:
        if record.id in antibody_ids:
            continue
        # left gaps, alignment, right gaps.
        match = re.match("(^-*)([^-].*[^-])(--*)$", str(record.seq))
        right_ends.append(match.start(3))
    maskpos = min(right_ends)

    # Cut the alignment down to just the antibodies and mask to just the
    # rightward region.
    antibodies_masked = NamedTemporaryFile("wt", buffering=1)
    for record in alignment:
        if record.id in antibody_ids:
            record_masked = SeqRecord(
                Seq("-" * maskpos + str(record.seq)[maskpos:]),
                id=record.id,
                description="")
            SeqIO.write(record_masked, antibodies_masked, "fasta")

    # Align the new alleles to this modified version.
    aligned_masked = NamedTemporaryFile("rt", buffering=1)
    shell(
        "clustalw -align -profile1={profile1} "
        "-profile2={profile2} -sequences -output=fasta "
        "-outfile={outfile}".format(
            profile1=antibodies_masked.name,
            profile2=alleles_fasta,
            outfile=aligned_masked.name))

    # Swap the original antibodies back into place.
    with open(output_fasta, "wt") as f_out:
        for record in alignment:
            if record.id in antibody_ids:
                SeqIO.write(record, f_out, "fasta")
        for record in SeqIO.parse(aligned_masked, "fasta"):
            if record.id in antibody_ids:
                continue
            SeqIO.write(record, f_out, "fasta")

def combine_aligned_segments(antibodies_fasta, with_v, with_d, with_j, output_fasta):
    """Combine the antibody sequences and separately aligned V(D)J alleles into one alignment."""
    antibodies = list(SeqIO.parse(antibodies_fasta, "fasta"))
    aligned_v = list(SeqIO.parse(with_v, "fasta"))
    aligned_d = list(SeqIO.parse(with_d, "fasta"))
    aligned_j = list(SeqIO.parse(with_j, "fasta"))
    antibody_ids = {record.id for record in antibodies}
    with open(output_fasta, "wt") as f_out:
        for record in aligned_v:
            SeqIO.write(record, f_out, "fasta")
        for record in aligned_d:
            if record.id not in antibody_ids:
                SeqIO.write(record, f_out, "fasta")
        for record in aligned_j:
            if record.id not in antibody_ids:
                SeqIO.write(record, f_out, "fasta")
