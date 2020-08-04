"""
Summarizing and reporting helper functions - IgDiscover.
"""

import csv
import re
import logging
from tempfile import NamedTemporaryFile
from snakemake.shell import shell
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from igseq.util import strdist_iupac_squeezed

LOGGER = logging.getLogger(__name__)

def align_next_segment(target_fasta, aligned_fasta, alleles_fasta, output_fasta):
    """Align discovered alleles from a segment to the regions not yet aligned to.

    target_fasta: path to FASTA for alignment target, i.e., known antibody
                  sequences
    aligned_fasta: path to FASTA for previous segment's alignment
    alleles_fasta: path to FASTA for discovered alleleles for next segment
    output_fasta: path to FASTA to save new alignment to
    """
    targets = list(SeqIO.parse(target_fasta, "fasta"))
    alignment = list(SeqIO.parse(aligned_fasta, "fasta"))
    target_ids = {record.id for record in targets}

    # Handle the case that we have no targets on record for this subject
    if not target_ids:
        shell("touch {output_fasta}")
        return

    # Find the farthest-left right edge of the previously-aligned set of alleles.
    right_ends = []
    for record in alignment:
        if record.id in target_ids:
            continue
        # left gaps, alignment, right gaps.
        match = re.match("(^-*)([^-].*[^-])(--*)$", str(record.seq))
        right_ends.append(match.start(3))
    maskpos = min(right_ends)

    # Cut the alignment down to just the targets and mask to just the
    # rightward region.
    targets_masked = NamedTemporaryFile("wt", buffering=1)
    for record in alignment:
        if record.id in target_ids:
            record_masked = SeqRecord(
                Seq("-" * maskpos + str(record.seq)[maskpos:]),
                id=record.id,
                description="")
            SeqIO.write(record_masked, targets_masked, "fasta")

    # Align the new alleles to this modified version.
    aligned_masked = NamedTemporaryFile("rt", buffering=1)
    shell(
        "clustalw -align -profile1={profile1} "
        "-profile2={profile2} -sequences -output=fasta "
        "-outfile={outfile}".format(
            profile1=targets_masked.name,
            profile2=alleles_fasta,
            outfile=aligned_masked.name))

    shell("cp {outfile_temp} {outfile_real}".format(
        outfile_temp=aligned_masked.name, outfile_real=output_fasta))

def combine_aligned_segments(target_fasta, with_v, with_d, with_j, output_fasta):
    """Combine the target sequences and separately aligned V(D)J alleles into one alignment."""
    targets = list(SeqIO.parse(target_fasta, "fasta"))
    if not targets:
        shell("touch {output_fasta}")
        return
    aligned = {
        "v": list(SeqIO.parse(with_v, "fasta")),
        "j": list(SeqIO.parse(with_j, "fasta"))
        }
    if with_d:
        aligned["d"] = list(SeqIO.parse(with_d, "fasta"))
    ab_ids = {record.id for record in targets}
    targets = {key: [entry for entry in aligned[key] if entry.id in ab_ids] for key in aligned}
    lengths = {key: max([len(rec) for rec in targets[key]]) for key in targets}
    if len(set(lengths)) > 1:
        LOGGER.warning("Aligned targets differ in length; will pad to max length.")
    with open(output_fasta, "wt") as f_out:
        for record in aligned["v"]:
            record.seq = Seq(str(record.seq).ljust(max(lengths.values()), "-"))
            SeqIO.write(record, f_out, "fasta")
        if "d" in aligned:
            for record in aligned["d"]:
                record.seq = Seq(str(record.seq).ljust(max(lengths.values()), "-"))
                if record.id not in ab_ids:
                    SeqIO.write(record, f_out, "fasta")
        for record in aligned["j"]:
            record.seq = Seq(str(record.seq).ljust(max(lengths.values()), "-"))
            if record.id not in ab_ids:
                SeqIO.write(record, f_out, "fasta")

def convert_combined_alignment(fasta_in, csv_out, specimens, antibody_isolates, wildcards):
    """Convert VDJ+antibody alignment from FASTA into CSV form."""
    fields = ["Specimen", "Subject", "Chain", "Type", "Category", "LineageDist", "SeqName", "Seq"]
    segments = ["IGHV", "IGHD", "IGHJ", "IGLV", "IGLJ", "IGKV", "IGKJ"]
    categories = ["???", "AntibodyLineage", "AntibodyIsolate"] + segments
    subject_lut = {specimen["Specimen"]: specimen["Subject"] for specimen in specimens.values()}
    rows = []
    def categorize(seqid):
        if seqid in antibody_isolates.keys():
            return "AntibodyIsolate"
        if seqid in {entry["AntibodyLineage"] for entry in antibody_isolates.values()}:
            return "AntibodyLineage"
        for segment in segments:
            match = re.match("^(" + segment + ").*$", seqid)
            if match:
                return match.group(1)
        return "???"
    for record in SeqIO.parse(fasta_in, "fasta"):
        row = {}
        row["Specimen"] = wildcards.specimen
        row["Subject"] = subject_lut[wildcards.specimen]
        row["Chain"] = wildcards.chain
        row["Type"] = wildcards.chain_type
        row["Category"] = categorize(record.id)
        row["SeqName"] = record.id
        row["Seq"] = str(record.seq)
        rows.append(row)
    lineage = ""
    for row in rows:
        if row["Category"] == "AntibodyLineage":
            if lineage:
                LOGGER.error("antibody lineage sequence already defined; check metadata.")
            lineage = row["Seq"]
    if not lineage:
        LOGGER.warning("No lineage sequence present; defaulting to first non-segment sequence.")
        for row in rows:
            if row["Category"] not in segments:
                lineage = row["Seq"]
    for row in rows:
        row["LineageDist"] = strdist_iupac_squeezed(row["Seq"], lineage)
    _write_alignment_csv(csv_out, rows, fields, categories)

def _write_alignment_csv(csv_out, rows, fields, categories):
    """Write alignment dictionaries as CSV, sorting by field/category."""
    # make a (sortable) list for each row dictionary by using the keys in their
    # given order.  Special case: category will be sorted according to an
    # explicitly defined order (R factor-like, kinda).
    def sortable_row(row):
        entries = []
        for field in fields:
            entry = row[field]
            if field == "Category":
                entry = categories.index(entry)
            entries.append(entry)
        return entries
    # Sort by the given ordered fields above
    rows = sorted(rows, key=sortable_row)
    with open(csv_out, "wt") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

def gather_antibodies(lineage, chain, antibody_isolates, output_fasta):
    """Gather antibody sequences from metadata for use in alignments."""
    LOGGER.info("gather_antibodies: lineage: %s", lineage)
    LOGGER.info("gather_antibodies: chain: %s", chain)
    LOGGER.info("gather_antibodies: antibody_isolates: %d", len(antibody_isolates))
    LOGGER.info("gather_antibodies: output_fasta: %s", output_fasta)
    is_lineage = lambda val: val["AntibodyLineageAttrs"]["AntibodyLineage"] == lineage
    isolates = {k: val for k, val in antibody_isolates.items() if is_lineage(val)}
    LOGGER.info("gather_antibodies: kept %d isolates matching lineage", len(isolates))
    lineages = {val["AntibodyLineage"]: val["AntibodyLineageAttrs"] for val in isolates.values()}
    LOGGER.info("gather_antibodies: kept %d entries matching lineage", len(lineages))
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
