"""
Helpers for running SONAR.
"""

import re
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from snakemake.io import expand
from . import data

LOGGER = logging.getLogger(__name__)

def gather_germline(input_fps, output_fp):
    """
    Prepare germline references for SONAR.

    Note that it expects a particular sequence naming scheme for genes and
    module 2 will fail if we don't follow it.
    """
    with open(output_fp, "w") as f_out:
        seqids_seen = set()
        for input_fp in input_fps:
            for record in SeqIO.parse(input_fp, "fasta"):
                new_id = munge_seqid_for_sonar(record.id, seqids_seen)
                record.id = new_id
                record.description = ""
                SeqIO.write(record, f_out, "fasta")

def munge_seqid_for_sonar(seqid, seqids_seen):
    """Modify given sequence ID to work with SONAR.

    SONAR expects a particular format for the sequence ID entries for V(D)J
    genes and alleles.  This will try to parse out  the right portion and
    return it.  If it can't be found, an error is logged and the original ID is
    returned instead.  This works with IMGT IDs and at least one other
    reference database we're using.

    Also, the given set of already-observed sequence IDs will be used to record
    each new ID and check against the previous set.  An error is logged if a
    sequence ID has already been observed.  Note that the set is modified in place.
    """
    # See lineage/2.1-calculate_id-div.py
    # (IG[HKL]V[^*]+)[^,\s]+
    match = re.search(r"\|?(IG[HKL][VDJ][-*0-9A-Za-z_]+)\|?", seqid)
    if match:
        new_id = match.group(1)
    else:
        LOGGER.error("Sequence ID format not recognized, SONAR might get upset: %s", seqid)
        new_id = seqid
    if new_id in seqids_seen:
        LOGGER.error("Sequence ID already seen in this set: %s", new_id)
    seqids_seen.add(new_id)
    return new_id

def gather_mature(antibody_isolates, subject, chain, output_fp):
    """Get the heavy or light antibody sequences for a given subject as FASTA."""
    col_key = {"heavy": "HeavySeq", "light": "LightSeq"}.get(chain)
    if not col_key:
        raise data.MetadataError("chain should be \"heavy\" or \"light\", not \"%s\"" % chain)
    for isolate_name, isolate_attrs in antibody_isolates.items():
        if isolate_attrs["Subject"] == subject:
            seq = isolate_attrs[col_key]
            seqid = isolate_name + "_" + col_key
            rec = SeqRecord(Seq(seq), id=seqid, description="")

def get_antibody_alleles(antibody_lineages, subject, chain, segment):
    """Get a list of antibody lineage germline allele IDs for a given subject/chain/segment.

    This could be more than one, but right now we're assuming elsewhere that
    we're only pursuing one lineage per subject.  A MetadataError is raised if
    an allele definition is blank.
    """
    seqids = []
    for lineage_attrs in antibody_lineages.values():
        if lineage_attrs["Subject"] == subject:
            key = segment + {"heavy": "H", "light": "L"}[chain]
            seqid = lineage_attrs[key]
            if not seqid:
                msg = "Missing allele definition for antibody lineage %s %s chain" % (
                    lineage_attrs["AntibodyLineage"], chain)
                raise data.MetadataError(msg)
            seqids.append(seqid)
    return seqids

def igdiscover_final_db(samples, subject, chain, chain_type, segment):
    """Get an IgDiscover output DB FASTA for a single segment."""
    if chain_type == "gamma":
        chain_type_naive = "mu"
    elif chain_type == "lambda":
        chain_type_naive = "lambda"
    elif chain_type == "kappa":
        chain_type_naive = "kappa"
    else:
        raise ValueError(
            "chain type should be gamma or lambda or kappa, not \"%s\"" % chain_type)
    # Again, here we want to refer to the naive specimens, not the IgG+ SONAR ones.
    specimen_match = lambda x: \
        x["Chain"] == chain and \
        x["Type"] == chain_type_naive and \
        "IgM+" in x["SpecimenAttrs"]["CellType"] and \
        x["SpecimenAttrs"]["Subject"] == subject
    specimens = [x["Specimen"] for x in samples.values() if specimen_match(x)]
    return expand(
        "analysis/igdiscover/{chain}.{chain_type}/{specimen}/final/database/{segment}.fasta",
        chain=chain,
        chain_type=chain_type_naive,
        specimen=specimens,
        segment=segment)

def sonar_module_1_inputs(wildcards):
    """List the inputs needed for SONAR module 1.

    This will be the V(D)J germline FASTA files and the FASTQ input.
    """
    if wildcards.chain == "heavy":
        segments = ["V", "D", "J"]
    elif wildcards.chain == "light":
        segments = ["V", "J"]
    else:
        raise ValueError('Chain should be either "heavy" or "light"')
    targets = dict(zip(segments, expand(
        "analysis/sonar/{subject}/{chain}.{chain_type}/germline.{segment}.fasta",
        subject=wildcards.subject,
        chain=wildcards.chain,
        chain_type=wildcards.chain_type,
        segment=segments)))
    targets["fastq"] = expand(
        "analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/{specimen}.fastq",
        subject=wildcards.subject,
        chain=wildcards.chain,
        chain_type=wildcards.chain_type,
        specimen=wildcards.specimen)[0]
    return targets

def sonar_module_2_inputs(wildcards):
    """List the inputs needed for SONAR module 2 (automatic)."""
    if wildcards.chain == "heavy":
        segments = ["V", "D", "J"]
    elif wildcards.chain == "light":
        segments = ["V", "J"]
    else:
        raise ValueError('Chain should be either "heavy" or "light"')
    targets = dict(zip(segments, expand(
        "analysis/sonar/{subject}/{chain}.{chain_type}/germline.{segment}.fasta",
        subject=wildcards.subject,
        chain=wildcards.chain,
        chain_type=wildcards.chain_type,
        segment=segments)))
    targets["module1"] = expand(
        "analysissonar/{subject}/{chain}.{chain_type}/{specimen}/"
        "output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        subject=wildcards.subject,
        chain=wildcards.chain,
        chain_type=wildcards.chain_type,
        specimen=wildcards.specimen)[0]
    targets["mab"] = expand(
        "analysis/sonar/{subject}/{chain}.{chain_type}/mab.fasta",
        subject=wildcards.subject,
        chain=wildcards.chain,
        chain_type=wildcards.chain_type,
        specimen=wildcards.specimen)[0]
    return targets

def sonar_islands_for_subject(samples, subject, chain, chain_type):
    """Get the islandSeqs FASTA files for all specimens for a particular subject/chain/chain_type"""
    islight = lambda s: "IgG+" in s["SpecimenAttrs"]["CellType"] and s["Type"] == chain_type
    specimens = [entry["Specimen"] for entry in samples.values() if islight(entry)]
    targets = {}
    for spec in specimens:
        targets[spec] = (
            "analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/"
            "output/sequences/nucleotide/{specimen}_islandSeqs.fa").format(
                subject=subject,
                chain=chain,
                chain_type=chain_type,
                specimen=specimens)
    return targets
