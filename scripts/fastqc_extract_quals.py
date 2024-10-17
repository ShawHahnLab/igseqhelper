#!/usr/bin/env python

"""
Extract per-base sequence quality data from FastQC .zip file as CSV.
"""

import re
import sys
import csv
import zipfile

def parse_data_for(handle, section):
    """Extract tabular data for named section from fastqc_data.txt lines"""
    chunk = []
    in_section = False
    for line in handle:
        line = line.decode("ASCII").rstrip("\n")
        if in_section:
            if line == ">>END_MODULE":
                break
            chunk.append(line.split("\t"))
        else:
            match = re.match(r">>(.*)\t(.*)", line)
            if match and match.group(1) == section:
                in_section = True
    return chunk

def _getname(zip_in, filename):
    for name in zip_in.namelist():
        if name.endswith(f"/{filename}"):
            break
    else:
        raise ValueError(f"{filename} not found")
    return name

def extract_fastqc_quals(path_zip, path_out):
    """Extract per-base sequence quality data from FastQC .zip file as CSV"""
    with zipfile.ZipFile(path_zip) as zip_in:
        with zip_in.open(_getname(zip_in, "fastqc_data.txt")) as f_in:
            data = parse_data_for(f_in, "Per base sequence quality")
        with open(path_out, "w", encoding="ASCII") as f_out:
            writer = csv.writer(f_out, lineterminator="\n")
            writer.writerows(data)

if __name__ == "__main__":
    extract_fastqc_quals(sys.argv[1], sys.argv[2])
