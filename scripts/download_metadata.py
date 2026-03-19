#!/usr/bin/env python

"""Download CSV(s) from google sheets.

YAML like this:
      google_slug: 2PACX-1vRsc...
      resources:
        - name: subjects
          path: subjects.csv
          google_gid: 1234567890
        - name: specimens
          path: specimens.csv
          google_gid: 9876543210
        ...
Becomes:
     metadata/samples.csv
     metadata/specimens.csv
     metadata/...

Based on the CSV URLs from Google that look like:

    https://docs.google.com/spreadsheets/d/e/{slug}/pub?gid={gid}&single=true&output=csv
"""

import csv
import argparse
import urllib.request
from io import StringIO
from pathlib import Path
import yaml

def download_metadata(path_yaml, path_out=None):
    """Download CSV from google sheets."""
    with open(path_yaml, encoding="ASCII") as f_in:
        info = yaml.safe_load(f_in)
    slug = info["google_slug"]
    for res in info["resources"]:
        gid = res["google_gid"]
        url = "https://docs.google.com/spreadsheets/d/e/" \
            f"{slug}/pub?gid={gid}&single=true&output=csv"
        parent = Path(path_yaml).parent
        path = parent / res["path"]
        if path_out is None or str(path) == path_out:
            with urllib.request.urlopen(url) as f_in, open(path, "wt", encoding="ASCII") as f_out:
                reader = csv.reader(StringIO(f_in.read().decode("UTF8")))
                writer = csv.writer(f_out, lineterminator="\n")
                for row in reader:
                    writer.writerow(row)

def main():
    """CLI for download_metadata"""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="YAML input")
    parser.add_argument("output", nargs="?", help="Specific CSV to download (default: all)")
    args = parser.parse_args()
    download_metadata(args.input, args.output)

if __name__ == "__main__":
    main()
