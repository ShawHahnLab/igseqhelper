"""
Some BLAST helpers.
"""

import csv
import gzip
from pathlib import Path
import logging
from functools import lru_cache
from snakemake.shell import shell
from Bio import SeqIO
from .taxonomy import get_lineage

LOGGER = logging.getLogger(__name__)

BLAST_COLS = "qseqid sseqid pident length mismatch gapopen qstart " \
    "qend sstart send evalue bitscore stitle sgi staxids"
BLAST_NT_DB = "/data/external/ncbi/blast/nt/nt"
BLAST_TAX_GROUPS = {
        40674: "Mammalia",
        9443: "Primates",
        9544: "Macaca mulatta",
        7962: "common carp", # probably just stubs that need trimming
        198431: "uncultured prokaryote",
        1922226: "Anderseniella sp. Alg231-50", # likely PhiX
        562: "Escherichia coli", # likely PhiX
        1224: "Proteobacteria",
        10239: "Viruses",
        11102: "HCV",
        10847: "PhiX"}

def blast(tsv_out, fasta_in, threads=1, dbpath=BLAST_NT_DB, blast_cols=BLAST_COLS):
    """Wrapper for blastn with TSV output."""
    cmd = "blastn -num_threads {threads} -query {fasta} -db '{db}' -outfmt '6 {cols}' -out '{out}'"
    cmd = cmd.format(
        threads=threads, fasta=fasta_in, db=dbpath, tsv_out=tsv_out, cols=blast_cols, out=tsv_out)
    shell(cmd)

def group_blast_results(csv_out, tsv_in, taxonomy_groups, taxids_skip=None, blast_cols=BLAST_COLS):
    """Assign the results in a BLAST TSV file based on the taxonomy IDs.

    This writes a CSV file with one or more rows for each row of the BLAST TSV
    file, assigning each staxid to a group.  taxonomy_groups should be a
    dictionary mapping tax IDs to group names.  The most specific rank matching
    will be used for the group, so the dictionary can contain entries from the
    same lineages (e.g. mammals and primates).  taxids_skip is a list of
    taxonomic IDs to skip.
    """
    LOGGER.info("group_blast_results: csv_out: %s", csv_out)
    LOGGER.info("group_blast_results: tsv_in: %s", tsv_in)
    LOGGER.info("group_blast_results: taxonomy groups: %d groups", len(taxonomy_groups))
    LOGGER.info("group_blast_results: taxids_skip: %s", taxids_skip)
    LOGGER.info("group_blast_results: blast_cols: %s", blast_cols)
    try:
        # is string?
        blast_cols = blast_cols.split()
    except AttributeError:
        # already list
        pass
    if not taxids_skip:
        taxids_skip = []
    # force the grouping keys to strings to match the lineage info
    taxonomy_groups = {str(key): val for key, val in taxonomy_groups.items()}
    # This helps speed things up by avoiding file I/O on every single entry.
    # The grand total of our lineages cache on disk is on a few MB so loading
    # any number should be safe.
    lineage = (lru_cache(maxsize=None))(lambda staxid: get_lineage(staxid)[::-1])
    qseqid_staxids = []
    qseqid = None
    with open(tsv_in) as f_in, open(csv_out, "w", buffering=1) as f_out:
        reader = csv.DictReader(f_in, delimiter="\t", fieldnames=blast_cols)
        writer = csv.DictWriter(f_out, fieldnames=["qseqid", "staxid", "group"])
        writer.writeheader()
        for row in reader:
            # go through each taxonomic level, from lowest to highest, and
            # stop on the first (most-specific) match
            # note that we may have more than on taxId per entry
            staxids = row["staxids"].split(";")
            for staxid in staxids:
                # Continuing an existing chunk for one sequence ID.  If we've
                # already handled this taxid, skip it.
                if qseqid == row["qseqid"]:
                    if staxid in qseqid_staxids:
                        continue
                    qseqid_staxids.append(staxid)
                # Starting a new chunk for a new sequence ID
                else:
                    qseqid = row["qseqid"]
                    qseqid_staxids = [staxid]
                if staxid in taxids_skip:
                    continue
                if staxid == "N/A":
                    group = "other"
                else:
                    for lvl in lineage(staxid):
                        group = taxonomy_groups.get(lvl["TaxId"])
                        if group:
                            break
                    else:
                        group = "other"
                out_row = {"qseqid": row["qseqid"], "staxid": staxid, "group": group}
                writer.writerow(out_row)

def load_blast_groupings(csv_in):
    """Load CSV from group_blast_results into a limited dictionary form."""
    with open(csv_in) as f_in:
        reader = csv.DictReader(f_in)
        last = None
        data = {}
        for row in reader:
            if last == row["qseqid"]:
                continue
            entry = {key: row[key] for key in row}
            del entry["qseqid"]
            data[row["qseqid"]] = entry
            last = row["qseqid"]
    return data

def group_seqs_by_blast(fqgz_in, dir_out, csv_in):
    """Divvy up sequences based on taxonomy of BLAST results.

    fqgz_in: path to one fastq.gz file for input
    dir_out: path to output directory.  individual per-group fastq.gz files
             will be written here.
    csv_in: path to CSV file created from group_blast_results for the input
            fastq.gz file.
    """
    LOGGER.info("group_seqs_by_blast: fqgz_in: %s", fqgz_in)
    LOGGER.info("group_seqs_by_blast: dir_out: %s", dir_out)
    LOGGER.info("group_seqs_by_blast: csv_in: %s", csv_in)
    groupings = load_blast_groupings(csv_in)
    LOGGER.info("group_seqs_by_blast: loaded %d seq assignments", len(groupings))
    writers = {}
    try:
        with gzip.open(fqgz_in, "rt") as f_in:
            for record in SeqIO.parse(f_in, "fastq"):
                try:
                    group = groupings[record.id]["group"]
                    LOGGER.debug("group_seqs_by_blast: record %s -> group %s", record.id, group)
                except KeyError:
                    group = "none"
                    LOGGER.debug("group_seqs_by_blast: record %s -> no group", record.id)
                try:
                    writer = writers[group]
                except KeyError:
                    fp_out = Path(dir_out) / (group + ".fastq.gz")
                    fp_out.parent.mkdir(parents=True, exist_ok=True)
                    LOGGER.debug("group_seqs_by_blast: setting up output: %s", fp_out)
                    writer = gzip.open(fp_out, "wt")
                    writers[group] = writer
                SeqIO.write(record, writer, "fastq")
    finally:
        for key in writers:
            writers[key].close()
