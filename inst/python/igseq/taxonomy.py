"""
Helpers for managing taxonomic IDs and lineages.
"""

import sys
import csv
import time
from subprocess import run
from threading import RLock
from pathlib import Path
from urllib.error import HTTPError
import logging
from Bio import Entrez

LINEAGE_FIELDNAMES = ["TaxId", "ScientificName", "Rank"]
LOGGER = logging.getLogger(__name__)
LINEAGE_LOCK = RLock()

# http://biopython.org/DIST/docs/tutorial/Tutorial.html
def download_lineage(taxid):
    """Downoad a list of the full lineage for the given taxid, including itself.

    The list is ordered from highest to lowest (e.g. kingdom to species)
    """
    if not Entrez.email:
        Entrez.email = run(["git", "config", "user.email"], check=True, capture_output=True).stdout
    tries = 3
    delay = 5
    while tries:
        try:
            handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
        except HTTPError as exception:
            # too many requests
            if exception.code == 429:
                LOGGER.warning(
                    "HTTP rate limit hit downloading lineage for %s; delaying %d s", taxid, delay)
                time.sleep(delay)
                tries -= 1
            else:
                LOGGER.critical("HTTP error downloading lineage for %s", taxid)
                raise exception
        else:
            break
    lineage = records[0]["LineageEx"]
    this = {key: records[0][key] for key in LINEAGE_FIELDNAMES}
    lineage.append(this)
    return lineage

def save_lineage(lineage, stream=sys.stdout):
    """Write a CSV file with the given lineage list."""
    writer = csv.DictWriter(stream, fieldnames=LINEAGE_FIELDNAMES)
    writer.writeheader()
    for entry in lineage:
        writer.writerow(entry)

def load_lineage(csv_in):
    """Read a CSV file into a lineage list."""
    with open(csv_in) as f_in:
        reader = csv.DictReader(f_in, fieldnames=LINEAGE_FIELDNAMES)
        return list(reader)

def get_lineage(taxid, cachedir="lineages"):
    """Get a list of the full lineage for the given taxid, including itself.

    This version caches requests as CSV files under a cache directry.
    """
    fp_csv = Path(cachedir) / (taxid + ".csv")
    if fp_csv.exists():
        return load_lineage(fp_csv)
    # If we need to update the cache, make sure only one thread at a time is
    # doing it!
    with LINEAGE_LOCK:
        lineage = download_lineage(taxid)
        fp_csv.parent.mkdir(parents=True, exist_ok=True)
        with open(fp_csv, "w") as f_out:
            save_lineage(lineage, f_out)
    return lineage
