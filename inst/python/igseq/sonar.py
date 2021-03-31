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

def setup_sonar_combos(sample_md_igg, antibody_lineages):
    """Create dictionary of lists of attributes used for various SONAR rules.

    The output is structured to work easily with expand(pattern, zip, **output)
    in snakemake rules.
    """
    lineages_per_subject = {}
    for lineage_name, lineage_attrs in antibody_lineages.items():
        if not lineage_attrs["Subject"] in lineages_per_subject:
            lineages_per_subject[lineage_attrs["Subject"]] = []
        lineages_per_subject[lineage_attrs["Subject"]].append(lineage_name)

    # These are the key attributes we refer to for SONAR, named to match the
    # wildcards used in the rules.
    combos = {
        "subject": [],
        "antibody_lineage": [],
        "chain": [],
        "chain_type": [],
        "specimen": []}
    for subject in lineages_per_subject:
        filter_to_subject = lambda vec: [
            x for x, y in zip(vec, sample_md_igg["subjects"]) if y == subject]
        for lineage in lineages_per_subject[subject]:
            chains = filter_to_subject(sample_md_igg["chains"])
            chain_types = filter_to_subject(sample_md_igg["chaintypes"])
            specimens = filter_to_subject(sample_md_igg["specimens"])
            for chain, chain_type, specimen in zip(chains, chain_types, specimens):
                combos["subject"].append(subject)
                combos["antibody_lineage"].append(lineage)
                combos["chain"].append(chain)
                combos["chain_type"].append(chain_type)
                combos["specimen"].append(specimen)
    return combos

def gather_germline(input_fps, output_fp):
    """
    Prepare germline references for SONAR.

    Note that it expects a particular sequence naming scheme for genes and
    module 2 will fail if we don't follow it.  Module 1 also seems confused by
    the periods IMGT uses in its sequences so we'll strip those out as well.
    """
    with open(output_fp, "w") as f_out:
        seqids_seen = set()
        for input_fp in input_fps:
            for record in SeqIO.parse(input_fp, "fasta"):
                # update with simple sequence ID
                new_id = munge_seqid_for_sonar(record.id, seqids_seen)
                record.id = new_id
                # update to remove periods in sequence content
                record.seq = Seq(str(record.seq).replace(".", ""))
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

def gather_mature(antibody_isolates, antibody_lineage, chain, output_fp):
    """Get the heavy or light antibody sequences for a given subject as FASTA."""
    col_key = {"heavy": "HeavySeq", "light": "LightSeq"}.get(chain)
    if not col_key:
        raise data.MetadataError("chain should be \"heavy\" or \"light\", not \"%s\"" % chain)
    with open(output_fp, "w") as f_out:
        for isolate_name, isolate_attrs in antibody_isolates.items():
            if isolate_attrs["AntibodyLineage"] == antibody_lineage:
                seq = isolate_attrs[col_key]
                seqid = isolate_name + "_" + col_key
                rec = SeqRecord(Seq(seq), id=seqid, description="")
                SeqIO.write(rec, f_out, "fasta")

def get_antibody_allele(antibody_lineages, antibody_lineage, subject, chain, segment):
    """Get the antibody lineage germline allele ID for a given lineage/subject/chain/segment.

    A MetadataError is raised if an allele definition is blank or if no
    matching lineage is found.
    """
    for lineage_attrs in antibody_lineages.values():
        if lineage_attrs["Subject"] == subject and \
                lineage_attrs["AntibodyLineage"] == antibody_lineage:
            key = segment + {"heavy": "H", "light": "L"}[chain]
            seqid = lineage_attrs[key]
            if not seqid:
                msg = "Missing allele definition for antibody lineage %s %s chain" % (
                    lineage_attrs["AntibodyLineage"], chain)
                raise data.MetadataError(msg)
            return seqid
    raise data.MetadataError("No antibody lineage %s found for subject %s" % (
        antibody_lineage, subject))

def igdiscover_final_db(samples, subject, chain, chain_type, segment):
    """Get an IgDiscover output DB FASTA for a single segment."""
    if chain == "heavy":
        chain_type_naive = "mu"
    elif chain == "light":
        if chain_type == "lambda":
            chain_type_naive = "lambda"
        elif chain_type == "kappa":
            chain_type_naive = "kappa"
        else:
            raise ValueError(
                "chain type should be gamma or lambda or kappa, not \"%s\"" % chain_type)
    else:
        raise ValueError("chain should be heavy or light, not \"%s\"" % chain)
    # Again, here we want to refer to the naive specimens, not the IgG+ SONAR ones.
    specimen_match = lambda x: \
        x["Chain"] == chain and \
        x["Type"] == chain_type_naive and \
        "IgM+" in x["SpecimenAttrs"]["CellType"] and \
        x["SpecimenAttrs"]["Subject"] == subject
    specimens = list({x["Specimen"] for x in samples.values() if specimen_match(x)})
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
    targets = expand(
        "analysis/sonar/{subject}/{chain}.{chain_type}/germline.{segment}.fasta",
        subject=wildcards.subject,
        chain=wildcards.chain,
        chain_type=wildcards.chain_type,
        segment=segments)
    targets = dict(zip(segments, targets))
    targets["fastq"] = expand(
        "analysis/sonar/{subject}/{chain}.{chain_type}/{antibody_lineage}/" \
        "{specimen}/{specimen}.fastq",
        subject=wildcards.subject,
        chain=wildcards.chain,
        chain_type=wildcards.chain_type,
        antibody_lineage=wildcards.antibody_lineage,
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
