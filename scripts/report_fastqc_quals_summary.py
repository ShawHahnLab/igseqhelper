#!/usr/bin/env python

"""
Make a per-run CSV summary of FastQC-supplied quality info.

This reads all available outputs from the fastqc_extract_quals rule and
produces a single summary CSV.  The quality summary values are the minimum
value (across bases) of the 25% percentile (across reads) of the observed
quality scores from the FASTQ data, as supplied by FastQC, rounded to the
nearest integer.  This way if the lower-quality reads at a particular cycle are
abnormally low, it should show up in the output value for that read file.

Output columns are:

 * Run: Run ID
 * R1Q: Quality summary for R1 reads, for the first half of the bases
 * I1Q: Quality summary for I1 reads, taken across all bases
 * R2Q: Quality summary for R2 reads, for the first half of the bases
"""

import re
import argparse
from csv import DictReader, writer
from pathlib import Path
from collections import defaultdict

ANALYSIS = Path("analysis")

def _parse_bases(txt):
    # sometimes it's an X-Y range instead of a single base, for the FASTQ files
    # with longer reads.  If not, this will just be the same number twice
    base1 = int(re.sub("-.*", "", txt))
    base2 = int(re.sub(".*-", "", txt))
    return (base1, base2)

def _min_val_in_range(quals, colname, min_base_fract=0, max_base_fract=1):
    num_bases = max(max(_parse_bases(row["#Base"])) for row in quals)
    min_base = 1 + min_base_fract*num_bases
    max_base = max_base_fract*num_bases
    vals = []
    for row in quals:
        base1, base2 = _parse_bases(row["#Base"])
        if min_base <= base2 and base1 <= max_base:
            vals.append(float(row[colname]))
    return min(vals)

def fastqc_quals_summary(path_out="/dev/stdout"):
    """Make a per-run summary of FastQC-supplied quality info."""
    quals = defaultdict(dict)
    for path_csv in (ANALYSIS/"fastqc/reads").glob("**/*fastqc.quals.csv"):
        runid = path_csv.parent.name
        rpkey = re.search(r"(..)_001_fastqc\.quals\.csv$", path_csv.name).group(1)
        with open(path_csv, encoding="ASCII") as f_in:
            if rpkey in quals[runid]:
                raise ValueError(f"Already have data for {runid} {rpkey}?")
            quals[runid][rpkey] = list(DictReader(f_in))
    out = []
    for runid, trio in quals.items():
        out.append([
            runid,
            round(_min_val_in_range(trio["R1"], "Lower Quartile", 0, 0.5)),
            round(_min_val_in_range(trio["I1"], "Lower Quartile")),
            round(_min_val_in_range(trio["R2"], "Lower Quartile", 0, 0.5))])
    out.sort()
    with open(path_out, "w", encoding="ASCII") as f_out:
        csvwriter = writer(f_out, lineterminator="\n")
        csvwriter.writerow(["Run", "R1Q", "I1Q", "R2Q"])
        csvwriter.writerows(out)

def main():
    """CLI for fastqc_quals_summary"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default="/dev/stdout", help="output CSV path")
    args = parser.parse_args()
    fastqc_quals_summary(args.output)

if __name__ == "__main__":
    main()
