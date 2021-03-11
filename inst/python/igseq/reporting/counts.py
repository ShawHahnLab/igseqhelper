"""
Summarizing and reporting helper functions - Sequence counts.
"""

import csv
import re
import logging
from collections import defaultdict
from igseq.data import load_csv, get_samples_per_run, MetadataError

LOGGER = logging.getLogger(__name__)

def counts_sample_summary(inputs, csv_out, samples):
    """Take a list of sequence counts for fastq.gz files and summarize per-sample."""

    rows = []
    for samp in samples.values():
        if not samp["Run"]:
            continue
        # Gather up the counts files for any one sample and sum the counts
        pattern = "analysis/counts/demux/{run}/[0-9]+/{sample}.I1.fastq.gz.counts".format(
            run=samp["Run"],
            sample=samp["Sample"])
        matched_paths = [path for path in inputs if re.match(pattern, path)]
        cts = 0
        for path in matched_paths:
            with open(path, "rt") as f_in:
                cts += int(f_in.read().strip())
        # For samples where we know how many cells went in, enter the number of
        # cells and the ratio of reads to cells.  If we don't know just leave
        # it blank.
        try:
            cellcount = samp["SpecimenAttrs"]["CellCount"]
            ratio = divide(cts, cellcount)
        except KeyError:
            cellcount = ""
            ratio = ""
        rows.append({
            "Run": samp["Run"],
            "Sample": samp["Sample"],
            "CellCount": cellcount,
            "NumSequences": cts,
            "Ratio": ratio})

    # Separately from the named samples, we have counts for how many read were
    # not assigned to samples, and how many of those mapped to the PhiX genome,
    # for each run associated with a sample.
    for run in {samp["Run"] for samp in samples.values() if samp["Run"]}:
        for sample in ["unassigned", "unassigned.phix"]:
            pattern = f"analysis/counts/demux/{run}/[0-9]+/{sample}.I1.fastq.gz.counts"
            matches = [path for path in inputs if re.match(pattern, path)]
            cts = 0
            for path in matches:
                with open(path, "rt") as f_in:
                    cts += int(f_in.read().strip())
            rows.append({
                "Run": run,
                "Sample": sample,
                "CellCount": "",
                "NumSequences": cts,
                "Ratio": ""})

    rows = sorted(rows, key=lambda row: (row["Run"], row["Sample"], row["NumSequences"]))
    with open(csv_out, "wt") as f_out:
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
        if row_in["Sample"] == "unassigned.phix":
            # this case is a subset of unassigned, handled below
            continue
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
        try:
            timepoint = str(int(specimens[specimen]["Timepoint"]))
        except ValueError:
            timepoint = ""
        counts.append({
            "Subject": specimens[specimen]["Subject"],
            "Timepoint": timepoint,
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
        r"analysis/counts/presto/data/([a-z]+)\.([a-z]+)/([A-Za-z0-9]+).R1.fastq.counts")
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
        r"analysis/counts/presto/assemble/([a-z]+)\.([a-z]+)/([A-Za-z0-9]+)_assemble-pass.fastq.counts")
    with open(output_fp, "wt") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(counts)

def counts_presto_qual_summary(input_fps, output_fp, specimens):
    """Take a list of per-specimen presto quality sequence count files and make a summary table."""
    counts, fieldnames = amplicon_summary(
        input_fps,
        specimens,
        r"analysis/counts/presto/qual/([a-z]+)\.([a-z]+)/([A-Za-z0-9]+)_quality-pass.fastq.counts")
    with open(output_fp, "wt") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(counts)

def counts_sonar_module1_summary(input_fps, output_fp, attrs):
    """Take a list of per-specimen SONAR rearrangment files and make a summary table."""
    fields_attrs = ["subject", "antibody_lineage", "chain", "chain_type", "specimen"]
    fields_cts = ["good", "nonproductive", "indel", "noCDR3", "noJ", "noV", "stop"]
    with open(output_fp, "wt", 1) as f_out:
        writer = csv.DictWriter(f_out, lineterminator="\n", fieldnames=fields_attrs + fields_cts)
        writer.writeheader()
        for idx, input_fp in enumerate(input_fps):
            attrs_here = {key: vec[idx] for key, vec in attrs.items()}
            cts = {key: 0 for key in fields_cts}
            with open(input_fp) as f_in:
                reader = csv.DictReader(f_in, delimiter="\t")
                for row in reader:
                    cts[row["status"]] += int(row["duplicate_count"])
            for key in cts:
                attrs_here[key] = cts[key]
            writer.writerow(attrs_here)

def divide(val1, val2, fmt="{:.6f}"):
    """Divide val1 by val2 as floats and return formatted string."""
    try:
        num = float(val1)/float(val2)
    except ValueError:
        num = ""
    else:
        num = fmt.format(num)
    return num
